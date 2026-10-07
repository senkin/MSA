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
