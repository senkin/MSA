"""CPU-only checks of the GPU (sample-transposed, batched) attribution path.

In GPU mode, run_NNLS.py replaces the per-sample greedy signature removal with one
batched NNLS call per greedy step across all samples. These tests pin that it gives the
same attributions as the per-sample CPU path:

* through the scipy fallback of batched_solve_and_score(), and
* through MSA's cuML wrapper (_GpuBatchedSolver), driven by the NumPy/SciPy stand-in
  for cupy/rmm/cuml in tests/fake_gpu.

Neither needs a GPU. The real cuML solver is checked by tests/gpu_acceptance.sh.
Run from the repository root:  python -m unittest discover -s tests -v
"""
import sys
import types
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
import run_NNLS as R  # noqa: E402

N_SAMPLES = 12
THRESHOLDS = [0.0, 1e-4, 1e-2]


def load_inputs():
    """65 SBS96 reference signatures and the first N_SAMPLES of the GPU test dataset."""
    sigs = pd.read_csv(REPO / 'signature_tables/sigProfiler_SBS_signatures.csv', index_col=[0, 1])
    muts = pd.read_csv(REPO / 'input_mutation_tables/SIM_test_GPU/WGS_SIM_test_GPU.96.csv', index_col=[0, 1])
    return sigs, muts.iloc[:, :N_SAMPLES].reindex(sigs.index)


def make_args(weak_threshold, use_gpu, gpu_batch_size=4096):
    return types.SimpleNamespace(optimise_signatures=True, optimisation_strategy='removal',
                                 weak_threshold=weak_threshold, strong_threshold=0.0,
                                 normalise_mutations=False, verbose=False,
                                 use_gpu=use_gpu, gpu_batch_size=gpu_batch_size)


def attribute(sigs, muts, args):
    weights, _, stat_info, _, _ = R.process_samples_batch(muts, sigs, list(range(sigs.shape[1])), args)
    return weights, stat_info


class AttributionTest(unittest.TestCase):

    def test_exact_mixture_is_recovered(self):
        sigs, _ = load_inputs()
        true_weights = pd.Series({'SBS1': 0.2, 'SBS5': 0.5, 'SBS13': 0.3})
        observed = (sigs[true_weights.index] @ true_weights * 2000).values
        with np.errstate(invalid='ignore'):  # exact fit: Jensen-Shannon takes sqrt of a rounding -0
            weights, mutations, fitted, residuals, stat_info = R.perform_signature_attribution(
                observed, sigs[true_weights.index])
        np.testing.assert_allclose(weights, true_weights.values, atol=1e-9)
        np.testing.assert_allclose(mutations, 2000 * true_weights.values, rtol=1e-9)
        np.testing.assert_allclose(fitted, observed, atol=1e-9)
        self.assertAlmostEqual(stat_info[0], observed.sum())  # mutational burden
        self.assertAlmostEqual(stat_info[4], 1.0)             # cosine similarity


class BatchedPathTest(unittest.TestCase):
    """Batched (GPU-path) greedy with the scipy fallback vs the per-sample CPU greedy."""

    @classmethod
    def setUpClass(cls):
        cls.sigs, cls.muts = load_inputs()
        cls.cpu = {t: attribute(cls.sigs, cls.muts, make_args(t, use_gpu=False)) for t in THRESHOLDS}

    def test_batched_matches_per_sample(self):
        # use_gpu with no GPU backend configured runs the batched path on scipy
        self.assertIsNone(R._GPU_SOLVER)
        for t in THRESHOLDS:
            with self.subTest(weak_threshold=t):
                weights, stat_info = attribute(self.sigs, self.muts, make_args(t, use_gpu=True))
                cpu_weights, cpu_stat_info = self.cpu[t]
                # At t=0 the two paths may keep different zero-weight signatures (the decision
                # is rounding noise), which leaves the fitted weights unchanged.
                np.testing.assert_allclose(weights.values, cpu_weights.values, rtol=0, atol=1e-9)
                np.testing.assert_allclose(stat_info.values, cpu_stat_info.values, rtol=1e-9, atol=1e-9)


class GpuWrapperTest(unittest.TestCase):
    """_GpuBatchedSolver (device residency, chunking, precision, provenance) on a fake GPU."""

    FAKE_GPU = str(REPO / 'tests' / 'fake_gpu')
    THRESHOLD = 1e-4

    @classmethod
    def setUpClass(cls):
        if any(m in sys.modules for m in ('cupy', 'rmm', 'cuml')):
            raise unittest.SkipTest('a GPU stack is already imported in this process')
        sys.path.insert(0, cls.FAKE_GPU)
        cls.sigs, cls.muts = load_inputs()
        cls.cpu_weights, _ = attribute(cls.sigs, cls.muts, make_args(cls.THRESHOLD, use_gpu=False))

    @classmethod
    def tearDownClass(cls):
        R._GPU_SOLVER = None
        R._GPU_INITIALISED = False
        R._GPU_STREAM = None
        R._GPU_POOL_BYTES = 0
        for name in [m for m in sys.modules if m.split('.')[0] in ('cupy', 'rmm', 'cuml')]:
            del sys.modules[name]
        sys.path.remove(cls.FAKE_GPU)

    def run_gpu(self, precision, gpu_batch_size):
        R.configure_gpu_backend(precision)
        try:
            weights, _ = attribute(self.sigs, self.muts,
                                   make_args(self.THRESHOLD, use_gpu=True, gpu_batch_size=gpu_batch_size))
            return weights, R._GPU_SOLVER
        finally:
            R._GPU_SOLVER = None

    def test_fp64_auto_batch_matches_cpu(self):
        weights, solver = self.run_gpu('fp64', gpu_batch_size=0)
        self.assertIsNotNone(solver.auto_batch_size)
        np.testing.assert_allclose(weights.values, self.cpu_weights.values, rtol=0, atol=1e-9)

    def test_fp64_small_batches_match_cpu(self):
        # 100 problems per call: every greedy step is split over many solver calls
        weights, _ = self.run_gpu('fp64', gpu_batch_size=100)
        np.testing.assert_allclose(weights.values, self.cpu_weights.values, rtol=0, atol=1e-9)

    def test_mixed_precision_close_to_cpu(self):
        # fp32 storage rounding may flip a decision that sits within ~1e-8 of the
        # threshold; allow at most one such sample here
        weights, _ = self.run_gpu('mixed', gpu_batch_size=0)
        differing = (np.abs(weights.values - self.cpu_weights.values) > 1e-3).any(axis=1)
        self.assertLessEqual(int(differing.sum()), 1)

    def test_provenance(self):
        _, solver = self.run_gpu('fp64', gpu_batch_size=0)
        info = solver.provenance()
        self.assertEqual(info['precision'], 'fp64')
        self.assertEqual(info['device']['name'], 'Fake GPU (NumPy)')
        self.assertEqual(info['versions']['cuml'], 'fake')
        self.assertIsNotNone(info['auto_batch_size'])
        self.assertIn('cuml_build', info)


if __name__ == '__main__':
    unittest.main()
