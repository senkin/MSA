"""The penalty-free removal trajectory reproduces the per-threshold greedy removal.

remove_weak_signatures_trajectory() records the full removal order once per sample;
select_from_trajectory() must then give, for any weak threshold, exactly the signatures
remove_weak_signatures() keeps at that threshold - checked for the default and test-profile
penalty ranges, down to the output tables. Run from the repository root:
    python -m unittest discover -s tests -v
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


if __name__ == '__main__':
    unittest.main()
