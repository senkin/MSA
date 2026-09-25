"""How attributions are scored and the optimal penalty is picked, and the run comparison
used by the GPU acceptance check. Run from the repository root:
    python -m unittest discover -s tests -v
"""
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
from common_methods import calculate_stat_scores  # noqa: E402
from calculate_optimal_penalties import calculate_optimal_penalty  # noqa: E402


class StatScoresTest(unittest.TestCase):

    def test_confusion_matrix_metrics(self):
        signatures = ['A', 'B', 'C']
        truth = pd.DataFrame([[0.5, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], columns=signatures)
        reco = pd.DataFrame([[0.7, 0.0, 0.3], [0.0, 1.0, 0.0], [0.2, 0.0, 0.8]], columns=signatures)
        # sample 1: A true positive, B false negative, C false positive
        # sample 2: B true positive, A and C true negatives
        # sample 3: C true positive, A false positive, B true negative -> TP=3, FN=1, FP=2, TN=3
        sensitivity, specificity, precision, accuracy, F1, MCC = calculate_stat_scores(signatures, reco, truth)
        self.assertAlmostEqual(sensitivity, 3 / 4)
        self.assertAlmostEqual(specificity, 3 / 5)
        self.assertAlmostEqual(precision, 3 / 5)
        self.assertAlmostEqual(accuracy, 6 / 9)
        self.assertAlmostEqual(F1, 6 / 9)
        self.assertAlmostEqual(MCC, (3 * 3 - 2 * 1) / math.sqrt(5 * 4 * 5 * 4))


class OptimalPenaltyTest(unittest.TestCase):
    """The rule that picks the penalty (default: specificity >= 0.95, then best sensitivity)."""

    WEAK = [0, 0.001, 0.01, 0.1]

    def tables(self, sensitivity, specificity):
        return (pd.DataFrame({0: sensitivity}, index=self.WEAK),
                pd.DataFrame({0: specificity}, index=self.WEAK))

    def test_best_sensitivity_among_specific_enough(self):
        penalty = calculate_optimal_penalty(*self.tables([0.9, 0.8, 0.7, 0.6], [0.5, 0.96, 0.99, 1.0]))
        self.assertEqual(penalty, (0.001, 0))

    def test_ties_go_to_the_smallest_penalty(self):
        penalty = calculate_optimal_penalty(*self.tables([0.9, 0.8, 0.8, 0.6], [0.5, 0.96, 0.99, 1.0]))
        self.assertEqual(penalty, (0.001, 0))

    def test_fallback_when_specificity_is_never_reached(self):
        # best sensitivity among the penalties with the highest specificity
        penalty = calculate_optimal_penalty(*self.tables([0.9, 0.8, 0.7, 0.6], [0.5, 0.6, 0.7, 0.7]))
        self.assertEqual(penalty, (0.01, 0))

    def test_prioritising_sensitivity(self):
        penalty = calculate_optimal_penalty(*self.tables([0.9, 0.8, 0.7, 0.6], [0.5, 0.96, 0.99, 1.0]),
                                            metric_to_prioritise='sensitivity', threshold=0.75)
        self.assertEqual(penalty, (0.001, 0))

    def test_undefined_sensitivity_gives_no_penalty(self):
        # a signature never present in the truth has no sensitivity: no penalty from it
        penalty = calculate_optimal_penalty(*self.tables([np.nan] * 4, [0.5, 0.96, 0.99, 1.0]))
        self.assertTrue(all(np.isnan(p) for p in penalty))


class CompareAttributionsTest(unittest.TestCase):
    """bin/compare_attributions.py, which gates tests/gpu_acceptance.sh."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='msa_compare_'))
        rng = np.random.default_rng(0)
        self.weights = pd.DataFrame(rng.dirichlet(np.ones(4), 10), columns=list('ABCD'),
                                    index=pd.Index(range(10), name='Sample'))
        for run in ('ref', 'cand'):
            (self.tmp / run / 'boot').mkdir(parents=True)
            self.write(run, self.weights)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, run, table, name='output_X_SBS_weights_table.csv'):
        table.to_csv(self.tmp / run / name)
        table.to_csv(self.tmp / run / 'boot' / 'output_X_SBS_1_weights_table.csv')

    def compare(self, *options):
        return subprocess.run([sys.executable, str(REPO / 'bin/compare_attributions.py'),
                               str(self.tmp / 'ref'), str(self.tmp / 'cand'), *options],
                              capture_output=True, text=True).returncode

    def test_identical_runs_pass_exact(self):
        self.assertEqual(self.compare('--exact'), 0)

    def test_differing_fraction_is_gated(self):
        changed = self.weights.copy()
        changed.iloc[3, 0] += 0.01  # 1 of 20 attributions (2 tables x 10 samples) differs
        self.write('cand', changed)
        self.assertEqual(self.compare('--max-differing-fraction', '0.1'), 0)
        self.assertEqual(self.compare('--max-differing-fraction', '0.02'), 1)
        self.assertEqual(self.compare('--exact'), 1)

    def test_missing_table_fails(self):
        (self.tmp / 'cand' / 'boot' / 'output_X_SBS_1_weights_table.csv').unlink()
        self.assertEqual(self.compare(), 1)


if __name__ == '__main__':
    unittest.main()
