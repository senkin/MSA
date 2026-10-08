"""Bootstrap summary tables (bin/make_bootstrap_tables.py) on a small hand-checkable fixture,
and the penalty selection of bin/calculate_optimal_penalties.py. Run from the repository root:
    python -m unittest discover -s tests -v
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / 'bin'
sys.path.insert(0, str(BIN))
from calculate_optimal_penalties import select_penalties  # noqa: E402
from make_bootstrap_tables import (calculate_confidence_intervals, calculate_signature_CPs,  # noqa: E402
                                   calculate_CPs_for_thresholds)
import run_NNLS as R  # noqa: E402

SIGNATURES = ['SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS5']  # sigTest SBS96 signatures
SAMPLES = ['s0', 's1', 's2']
N_BOOTSTRAP = 4


def weights(rows):
    return pd.DataFrame(rows, index=pd.Index(SAMPLES, name='Sample'), columns=SIGNATURES, dtype=float)


# Bootstrap weights: SBS1 in every iteration, SBS2 in half of them, SBS4 always in sample s0
def bootstrap_weights(i):
    sbs2 = 0.4 if i <= 2 else 0.0
    return weights([[0.5, sbs2, 0, 0.1, 0], [0.5, sbs2, 0, 0, 0], [0.5, sbs2, 0, 0, 0]])


CENTRAL = weights([[0.6, 0.3, 0, 0.1, 0], [0.7, 0.3, 0, 0, 0], [0.6, 0.4, 0, 0, 0]])
# Truth: SBS1-3 act, SBS4 and SBS5 do not
TRUTH = weights([[0.5, 0.5, 0, 0, 0], [1.0, 0, 0, 0, 0], [0.5, 0, 0.5, 0, 0]])


class BootstrapTablesTest(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='msa_tables_test_'))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write_attributions(self, dataset, truth=False):
        folder = self.tmp / 'in' / dataset
        (folder / 'bootstrap_output').mkdir(parents=True)

        def write(prefix, table):
            table.to_csv(f'{prefix}_weights_table.csv')
            (table * 1000).to_csv(f'{prefix}_mutations_table.csv')
            pd.DataFrame({'Cosine': 0.9}, index=table.index).to_csv(f'{prefix}_stat_info.csv')

        write(folder / f'output_{dataset}_SBS', CENTRAL)
        for i in range(1, N_BOOTSTRAP + 1):
            write(folder / f'bootstrap_output/output_{dataset}_SBS_{i}', bootstrap_weights(i))
        if truth:
            TRUTH.to_csv(folder / f'WGS_{dataset}.96.weights.csv')

    def make_tables(self, dataset):
        env = dict(os.environ, MPLCONFIGDIR=str(self.tmp / '.matplotlib'), MPLBACKEND='Agg')
        result = subprocess.run(
            [sys.executable, str(BIN / 'make_bootstrap_tables.py'), '-d', dataset, '-t', 'SBS', '-c', '96',
             '-p', 'sigTest', '-S', str(REPO / 'signature_tables'), '-i', str(self.tmp / 'in'),
             '-o', str(self.tmp / 'out'), '-n', str(N_BOOTSTRAP), '-T', '0', '5', '10'],
            cwd=self.tmp, env=env, capture_output=True, text=True)
        if result.returncode != 0:
            raise AssertionError(f'make_bootstrap_tables.py failed:\n{result.stdout}\n{result.stderr}')
        return self.tmp / 'out' / dataset

    def test_real_data_named_with_SIM(self):
        # no truth table: real data, whatever the name
        self.write_attributions('MYSIMDATA')
        out = self.make_tables('MYSIMDATA')
        self.assertFalse((out / 'truth_studies').exists())
        prevalences = pd.read_csv(out / 'signatures_prevalences_MYSIMDATA_SBS.csv', index_col=0)
        np.testing.assert_array_equal(prevalences['SBS1'], 1)
        np.testing.assert_array_equal(prevalences['SBS2'], 0.5)
        np.testing.assert_array_equal(prevalences['SBS4'], [1, 0, 0])
        # pruned: zero where the lower CI bound is zero, the rounded central count elsewhere
        pruned = pd.read_csv(out / 'pruned_attribution_MYSIMDATA_SBS_abs_mutations.csv', index_col=0)
        np.testing.assert_array_equal(pruned['SBS1'], [600, 700, 600])
        np.testing.assert_array_equal(pruned['SBS2'], 0)
        np.testing.assert_array_equal(pruned['SBS4'], [100, 0, 0])

    def test_simulated_data_found_by_truth_table(self):
        self.write_attributions('TRUTHDATA', truth=True)
        out = self.make_tables('TRUTHDATA')
        cps = pd.read_csv(out / 'truth_studies/signatures_CPs_SBS.csv', index_col=0)
        # SBS1: CI [0.5, 0.5] holds the truth in s0 and s2; SBS3 (never attributed) in s0 and s1
        self.assertAlmostEqual(cps.loc[0, 'SBS1'], 2 / 3)
        self.assertAlmostEqual(cps.loc[0, 'SBS3'], 2 / 3)
        # Others pools SBS4 and SBS5: s0 covers 1 of 2 (SBS4 attributed), s1 and s2 both
        self.assertAlmostEqual(cps.loc[0, 'Others'], (0.5 + 1 + 1) / 3)
        for name in ('stat_scores_tables_SBS.csv', 'sensitivity_thresholds_SBS.csv'):
            self.assertTrue((out / 'truth_studies' / name).exists(), name)
        # summed over the signatures: truth minus the bootstrap mean / median, and minus the
        # central attribution (e.g. s0: SBS2 0.5 - 0.2, SBS4 0 - 0.1)
        difference = pd.read_csv(out / 'truth_studies/truth_and_measured_difference_SBS.csv', index_col=0)
        np.testing.assert_allclose(difference['Truth - mean'], [0.2, 0.3, 0.3])
        np.testing.assert_allclose(difference['Truth - median'], [0.2, 0.3, 0.3])
        np.testing.assert_allclose(difference['Truth - NNLS'], 0, atol=1e-12)


class TruthStudiesTest(unittest.TestCase):
    """Confidence intervals and their coverage of the truth (make_bootstrap_tables.py)."""

    # 2 samples x 3 signatures
    LOWER = np.array([[0, 0.2, 0], [0, 0, 0.1]])
    UPPER = np.array([[0.01, 0.3, 0], [0.2, 0, 0.3]])
    TRUTH = np.array([[0.03, 0.5, 0], [0.1, 0.04, 0.2]])
    SIGNATURES = ['A', 'B', 'C']

    def test_confidence_intervals(self):
        # 1 sample, 2 signatures, 4 bootstrap iterations
        attributions = np.array([[[0.4, 0.4, 0, 0], [0.5, 0.5, 0.5, 0.5]]])
        central = pd.DataFrame([[0.3, 0.5]], index=['s0'], columns=['A', 'B'])
        table, lower, upper = calculate_confidence_intervals(attributions, central, ['s0'], ['A', 'B'], 0.95)
        np.testing.assert_allclose(lower, [[0, 0.5]])
        np.testing.assert_allclose(upper, [[0.4, 0.5]])
        central_value, (cell_lower, cell_upper) = table.loc['s0', 'A']
        self.assertEqual((central_value, cell_lower, cell_upper), (0.3, 0, 0.4))

    def test_CPs_for_thresholds(self):
        CPs = calculate_CPs_for_thresholds(self.LOWER, self.UPPER, self.TRUTH, ['s0', 's1'],
                                           self.SIGNATURES, [0, 0.05])
        # threshold 0: s0 A (0.03 above [0, 0.01]) and s1 B (0.04 above [0, 0]) are missed;
        # s0 B is covered as CI and truth are both above zero
        np.testing.assert_allclose(CPs[0].loc[0], [0.5, 0.5, 1])
        # a CI from zero reaches up to the threshold: 0.03 and 0.04 are within 0.05
        np.testing.assert_allclose(CPs[0.05].loc[0], [1, 1, 1])

    def test_signature_CPs(self):
        CPs = calculate_signature_CPs(self.LOWER, self.UPPER, self.TRUTH, ['s0', 's1'],
                                      self.SIGNATURES, acting_signatures=['A', 'B'])
        # A: covered in s1 only; B: never (no threshold here, 0.5 is outside [0.2, 0.3])
        self.assertAlmostEqual(CPs.loc[0, 'A'], 0.5)
        self.assertAlmostEqual(CPs.loc[0, 'B'], 0)
        # Others pools the signatures absent from the simulation: C, covered in both samples
        self.assertAlmostEqual(CPs.loc[0, 'Others'], 1)


class SelectPenaltiesTest(unittest.TestCase):
    """The penalty pair passed to the final attribution."""

    WEAK, STRONG = [0, 0.001, 0.01], [0, 0.05]

    def table(self, rows):
        return pd.DataFrame(rows, index=['SBS1', 'SBS2', 'SBS3'],
                            columns=['optimal_weak_penalty', 'optimal_strong_penalty'], dtype=float)

    def select(self, table, average=(0.001, 0), **options):
        return select_penalties(table, average, self.WEAK, self.STRONG, **options)

    def test_most_conservative_signature(self):
        table = self.table([[0.001, 0], [0.01, 0], [np.nan, np.nan]])
        self.assertEqual(self.select(table), (0.01, 0))

    def test_prioritised_and_average(self):
        table = self.table([[0.001, 0], [0.01, 0], [0, 0.05]])
        self.assertEqual(self.select(table, signatures_to_prioritise=['SBS1']), (0.001, 0))
        self.assertEqual(self.select(table, average=(0, 0.05), use_average=True), (0, 0.05))

    def test_undefined_penalty_falls_back_with_a_warning(self):
        # e.g. every signature deprioritised: never pass NaN on to the final attribution
        table = self.table([[np.nan, np.nan]] * 3)
        with self.assertWarnsRegex(UserWarning, 'smallest scanned'):
            self.assertEqual(self.select(table), (0, 0))
        table = self.table([[0.01, 0], [np.nan, np.nan], [np.nan, np.nan]])
        with self.assertWarnsRegex(UserWarning, 'smallest scanned'):
            self.assertEqual(self.select(table, signatures_to_prioritise=['SBS2', 'SBS3']), (0, 0))
        with warnings.catch_warnings():
            warnings.simplefilter('error')  # nothing to warn about with a defined penalty
            self.assertEqual(self.select(table), (0.01, 0))


class ThresholdArgumentsTest(unittest.TestCase):
    """run_NNLS.py -W/-S values."""

    def test_parsed(self):
        self.assertEqual(R.parse_thresholds(['0.01', '0', '0.01']), [0.0, 0.01])
        self.assertEqual(R.parse_thresholds(['auto'], allow_auto=True), 'auto')

    def test_nan_rejected(self):
        # a NaN threshold would silently skip the optimisation
        with self.assertRaises(ValueError):
            R.parse_thresholds(['nan'])


if __name__ == '__main__':
    unittest.main()
