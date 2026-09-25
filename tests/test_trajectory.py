"""The penalty-free removal trajectory reproduces the per-threshold greedy removal.

remove_weak_signatures_trajectory() records the full removal order once per sample;
select_from_trajectory() must then give, for any weak threshold, exactly the signatures
remove_weak_signatures() keeps at that threshold - checked for the default and test-profile
penalty ranges, down to the output tables. The same holds for the batched (GPU-path)
trajectory against remove_weak_signatures_batched(), on the scipy fallback and on the
fake GPU (tests/fake_gpu), and for run_NNLS.py runs with several thresholds against one
run per threshold. Run from the repository root:
    python -m unittest discover -s tests -v
"""
import os
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
import run_NNLS as R  # noqa: E402

N_SAMPLES = 8
RANGES = {'default': [0, 0.0001, 0.0002, 0.0003], 'test profile': [0, 0.01, 0.02]}


class TrajectoryTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sigs = pd.read_csv(REPO / 'signature_tables/sigProfiler_SBS_signatures.csv', index_col=[0, 1])
        muts = pd.read_csv(REPO / 'input_mutation_tables/SIM_test_GPU/WGS_SIM_test_GPU.96.csv', index_col=[0, 1])
        cls.muts = muts.iloc[:, :N_SAMPLES].reindex(cls.sigs.index)
        cls.names = list(cls.sigs.columns)
        cls.sig_values = np.asarray(cls.sigs.values, dtype=float)
        cls.all_cols = list(range(len(cls.names)))
        cls.trajectories = {}
        for sample in cls.muts.columns:
            observed = np.asarray(cls.muts[sample].values, dtype=float)
            cls.trajectories[sample] = R.remove_weak_signatures_trajectory(
                observed, np.linalg.norm(observed), cls.sig_values, cls.all_cols)

    def test_matches_per_threshold_removal(self):
        for label, thresholds in RANGES.items():
            for w in thresholds:
                for sample in self.muts.columns:
                    with self.subTest(range=label, weak_threshold=w, sample=sample):
                        observed = np.asarray(self.muts[sample].values, dtype=float)
                        expected, _ = R.remove_weak_signatures(observed, np.linalg.norm(observed), self.sig_values,
                                                               self.all_cols, self.names, w)
                        active, trajectory = self.trajectories[sample]
                        self.assertEqual(R.select_from_trajectory(active, trajectory, w), expected)

    def test_output_tables_identical(self):
        def from_trajectory(selected_mutations, initial_signatures, all_available_signatures,
                            strategy='removal', weak_threshold=0.01, **kwargs):
            sample = next(s for s in self.muts.columns
                          if np.array_equal(self.muts[s].values, selected_mutations))
            keep = R.select_from_trajectory(*self.trajectories[sample], weak_threshold)
            return all_available_signatures.loc[:, [self.names[c] for c in keep]]

        args = types.SimpleNamespace(optimise_signatures=True, optimisation_strategy='removal',
                                     weak_threshold=0.0001, strong_threshold=0.0, normalise_mutations=False,
                                     verbose=False, use_gpu=False, gpu_batch_size=0)
        expected = R.process_samples_batch(self.muts, self.sigs, self.all_cols, args)
        original = R.optimise_signatures
        R.optimise_signatures = from_trajectory
        try:
            actual = R.process_samples_batch(self.muts, self.sigs, self.all_cols, args)
        finally:
            R.optimise_signatures = original
        for exp, act in zip(expected, actual):  # weights, mutations, stat_info, residuals, fitted
            self.assertTrue(np.array_equal(exp.values, act.values, equal_nan=True))

    def test_survivors_shrink_with_threshold(self):
        grid = [0.0] + list(np.logspace(-6, -1, 26))
        for sample, (active, trajectory) in self.trajectories.items():
            self.assertEqual(len(trajectory), len(active) - 1)
            kept = [set(R.select_from_trajectory(active, trajectory, w)) for w in grid]
            for looser, stricter in zip(kept, kept[1:]):
                self.assertLessEqual(stricter, looser, f'sample {sample}')


def gpu_args(weak_threshold, gpu_batch_size=0):
    return types.SimpleNamespace(optimise_signatures=True, optimisation_strategy='removal',
                                 weak_threshold=weak_threshold, strong_threshold=0.0,
                                 normalise_mutations=False, verbose=False,
                                 use_gpu=True, gpu_batch_size=gpu_batch_size, threshold_scan='trajectory')


class BatchedTrajectoryTest(unittest.TestCase):
    """GPU path: batched trajectory vs remove_weak_signatures_batched(), output tables bitwise."""

    FAKE_GPU = str(REPO / 'tests' / 'fake_gpu')
    THRESHOLDS = [0, 0.0001, 0.01]  # kept short for CI; wider coverage in the GPU acceptance check

    @classmethod
    def setUpClass(cls):
        cls.sigs = pd.read_csv(REPO / 'signature_tables/sigProfiler_SBS_signatures.csv', index_col=[0, 1])
        muts = pd.read_csv(REPO / 'input_mutation_tables/SIM_test_GPU/WGS_SIM_test_GPU.96.csv', index_col=[0, 1])
        cls.muts = muts.iloc[:, :6].reindex(cls.sigs.index)
        cls.muts.iloc[:, 2] = 0  # a sample without mutations is skipped by both paths
        cls.sel = list(range(cls.sigs.shape[1]))

    def assert_equivalent(self, gpu_batch_size):
        trajectories = R.removal_trajectories(self.muts, self.sigs, self.sel, gpu_args(None, gpu_batch_size))
        self.assertNotIn(2, trajectories)
        for w in self.THRESHOLDS:
            with self.subTest(weak_threshold=w, gpu_batch_size=gpu_batch_size):
                args = gpu_args(w, gpu_batch_size)
                expected = R.process_samples_batch(self.muts, self.sigs, self.sel, args)
                actual = R.process_samples_batch(self.muts, self.sigs, self.sel, args, trajectories)
                for exp, act in zip(expected, actual):
                    self.assertTrue(np.array_equal(exp.values, act.values, equal_nan=True))

    def test_scipy_fallback(self):
        self.assertIsNone(R._GPU_SOLVER)
        self.assert_equivalent(gpu_batch_size=4096)
        self.assert_equivalent(gpu_batch_size=100)  # every greedy step split over many solver calls

    def test_fake_gpu(self):
        if any(m in sys.modules for m in ('cupy', 'rmm', 'cuml')):
            self.skipTest('a GPU stack is already imported in this process')
        sys.path.insert(0, self.FAKE_GPU)
        try:
            for precision in ('fp64', 'mixed'):
                R.configure_gpu_backend(precision)
                try:
                    self.assert_equivalent(gpu_batch_size=0)
                finally:
                    R._GPU_SOLVER = None
        finally:
            R._GPU_INITIALISED, R._GPU_STREAM, R._GPU_POOL_BYTES = False, None, 0
            for name in [m for m in sys.modules if m.split('.')[0] in ('cupy', 'rmm', 'cuml')]:
                del sys.modules[name]
            sys.path.remove(self.FAKE_GPU)


class AutoThresholdsTest(unittest.TestCase):

    def test_threshold_labels(self):
        # labels name the output folders: existing ones must not change, new ones must be exact
        for value, label in [(0, '0'), (0.0, '0'), (0.0001, '0.0001'), (0.0003, '0.0003'),
                             (0.01, '0.01'), (0.02, '0.02'), (1, '1'), (1.2e-07, '1.2e-07')]:
            self.assertEqual(R.format_threshold(value), label)
        for value in (3.7e-05, 0.00012345, 5e-07, 0.18):
            self.assertEqual(float(R.format_threshold(value)), value)

    def test_breakpoint_based_selection(self):
        # contributions are not monotonic: the breakpoints are their running maxima
        trajectories = [[(0, 1e-17), (1, 2e-5), (2, 1e-5), (3, 3e-3)], [(4, -1e-16), (5, 0.2)], []]
        thresholds = R.auto_weak_thresholds(trajectories, n_thresholds=5)
        self.assertEqual(thresholds[0], 0.0)
        self.assertEqual(thresholds, sorted(set(thresholds)))
        self.assertLessEqual(len(thresholds), 6)
        self.assertGreaterEqual(min(thresholds[1:]), R.AUTO_THRESHOLD_FLOOR)  # noise ignored
        for w in thresholds:  # labels in output names map back to the same value
            self.assertEqual(float(R.format_threshold(w)), w)
        self.assertEqual(R.auto_weak_thresholds([[(0, 1e-18)]]), [0.0])


class MultiThresholdCliTest(unittest.TestCase):
    """run_NNLS.py with several thresholds writes exactly what one run per threshold writes."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='msa_multi_threshold_'))
        shutil.copytree(REPO / 'input_mutation_tables/SIM_test', self.tmp / 'inputs/SIM_test')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_nnls(self, *args):
        env = dict(os.environ, MPLCONFIGDIR=str(self.tmp / '.matplotlib'))
        result = subprocess.run([sys.executable, str(REPO / 'bin/run_NNLS.py'), '-d', 'SIM_test', '-t', 'SBS',
                                 '-c', '96', '-p', 'sigTest', '-i', str(self.tmp / 'inputs'),
                                 '-s', str(REPO / 'signature_tables'), '-x', '--add_suffix', *map(str, args)],
                                cwd=self.tmp, env=env, capture_output=True, text=True)
        if result.returncode != 0:
            raise AssertionError(f'run_NNLS.py failed:\n{result.stdout}\n{result.stderr}')

    def test_matches_one_run_per_threshold(self):
        thresholds = ['0', '0.01', '0.02']
        for w in thresholds:
            self.run_nnls('-W', w, '-S', 0, '-o', self.tmp / 'single')
            self.run_nnls('-W', w, '-S', 0, '-o', self.tmp / 'single', '-B', '--n_bootstrap', 2, '--seed', 3,
                          '--bootstrap_output_suffix', f'{w}_0')
        for scan in ('trajectory', 'separate'):
            with self.subTest(threshold_scan=scan):
                out = self.tmp / scan
                self.run_nnls('-W', *thresholds, '-S', 0, '-o', out, '--threshold_scan', scan)
                self.run_nnls('-W', *thresholds, '-S', 0, '-o', out, '--threshold_scan', scan,
                              '-B', '--n_bootstrap', 2, '--seed', 3)
                single = {p.relative_to(self.tmp / 'single'): p.read_bytes()
                          for p in sorted((self.tmp / 'single').rglob('*')) if p.is_file()}
                multi = {p.relative_to(out): p.read_bytes() for p in sorted(out.rglob('*')) if p.is_file()}
                self.assertEqual(len(single), 3 * (5 + 2 * 3))
                self.assertEqual(single, multi)

    def test_auto_resolves_labels(self):
        self.run_nnls('-W', 'auto', '-S', 0, '--auto_thresholds', 3, '--resolve_thresholds_only',
                      '--thresholds_dir', self.tmp / 'labels', '-o', self.tmp / 'unused')
        labels = (self.tmp / 'labels/weak_thresholds.txt').read_text().split()
        self.assertEqual(labels[0], '0')
        self.assertTrue(2 <= len(labels) <= 4)
        self.assertEqual((self.tmp / 'labels/strong_thresholds.txt').read_text().split(), ['0'])
        self.assertFalse((self.tmp / 'unused').exists())


if __name__ == '__main__':
    unittest.main()
