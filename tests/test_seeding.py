"""Reproducibility of the pipeline's random steps (simulations, bootstrap resampling).

Runs the scripts the way the Nextflow processes do. Run from the repository root:
    python -m unittest discover -s tests -v
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / 'bin'
SIGNATURES = REPO / 'signature_tables'
sys.path.insert(0, str(BIN))
import run_NNLS as R  # noqa: E402


def run_script(script, *args, cwd):
    env = dict(os.environ, MPLCONFIGDIR=str(Path(cwd) / '.matplotlib'), MPLBACKEND='Agg')
    result = subprocess.run([sys.executable, str(BIN / script), *map(str, args)],
                            cwd=cwd, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError(f'{script} failed:\n{result.stdout}\n{result.stderr}')


def read_bytes(folder):
    """{relative path: bytes} for every file under folder."""
    folder = Path(folder)
    return {str(p.relative_to(folder)): p.read_bytes() for p in sorted(folder.rglob('*')) if p.is_file()}


class SeedingTest(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='msa_seed_test_'))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def bootstrap(self, out, start, n, seed):
        # same arguments as the run_bootstrap_NNLS / run_final_bootstrap_NNLS processes
        run_script('run_NNLS.py', '-B', '--n_bootstrap', n, '--bootstrap_start_index', start,
                   '--bootstrap_method', 'binomial', '-d', 'SIM_test', '-t', 'SBS', '-c', 96,
                   '-x', '-W', 0.01, '-S', 0, '-p', 'sigTest', '-i', self.tmp / 'inputs',
                   '-s', SIGNATURES, '-o', out, '--seed', seed, cwd=self.tmp)

    def test_bootstrap_independent_of_chunking(self):
        shutil.copytree(REPO / 'input_mutation_tables/SIM_test', self.tmp / 'inputs/SIM_test')
        self.bootstrap(self.tmp / 'one_task', start=1, n=4, seed=11)
        self.bootstrap(self.tmp / 'two_tasks', start=1, n=2, seed=11)
        self.bootstrap(self.tmp / 'two_tasks', start=3, n=2, seed=11)
        one, two = read_bytes(self.tmp / 'one_task'), read_bytes(self.tmp / 'two_tasks')
        self.assertEqual(len(one), 12)  # 4 iterations x (weights, mutations, stat_info)
        self.assertEqual(one, two)

        self.bootstrap(self.tmp / 'other_seed', start=1, n=1, seed=12)
        name = 'SIM_test/bootstrap_output/output_SIM_test_SBS_1_mutations_table.csv'
        self.assertNotEqual(read_bytes(self.tmp / 'other_seed')[name], one[name])

    def simulate(self, out, seed):
        run_script('simulate_data.py', '-d', 'SIM_seed', '-t', 'SBS', '-c', 96, '-n', 20,
                   '-p', 'sigTest', '-s', SIGNATURES, '-z', '--noise_type', 'gaussian', '-Z', 10,
                   '-o', out, '--seed', seed, cwd=self.tmp)
        folder = Path(out) / 'SIM_seed'
        return [(folder / f).read_bytes() for f in ('WGS_SIM_seed.96.csv', 'WGS_SIM_seed.96.weights.csv')]

    def test_simulation_seeded(self):
        first = self.simulate(self.tmp / 'a', seed=7)
        self.assertEqual(first, self.simulate(self.tmp / 'b', seed=7))
        self.assertNotEqual(first, self.simulate(self.tmp / 'c', seed=8))

        counts = pd.read_csv(self.tmp / 'a/SIM_seed/WGS_SIM_seed.96.csv', index_col=[0, 1])
        weights = pd.read_csv(self.tmp / 'a/SIM_seed/WGS_SIM_seed.96.weights.csv', index_col=0)
        self.assertEqual(counts.shape[1], 20)
        self.assertEqual(len(weights), 20)
        self.assertTrue((counts.values >= 0).all() and np.issubdtype(counts.values.dtype, np.integer))
        row_sums = weights.sum(axis=1).values
        self.assertTrue(np.all(row_sums <= 1 + 1e-9))  # true weights, minus zero-inflated ones


class BootstrapMethodsTest(unittest.TestCase):

    def setUp(self):
        self.table = pd.read_csv(REPO / 'input_mutation_tables/SIM_test/WGS_SIM_test.96.csv', index_col=[0, 1])
        self.table.iloc[:, 1] = 0  # a sample without mutations
        np.random.seed(0)

    def test_multinomial_keeps_totals(self):
        resampled = R.bootstrap_mutation_table(self.table, 'multinomial')
        self.assertTrue(resampled.index.equals(self.table.index) and resampled.columns.equals(self.table.columns))
        np.testing.assert_array_equal(resampled.sum().values, self.table.sum().values)

    def test_binomial_gives_counts(self):
        resampled = R.bootstrap_mutation_table(self.table, 'binomial')
        self.assertTrue(resampled.index.equals(self.table.index) and resampled.columns.equals(self.table.columns))
        self.assertTrue((resampled.values >= 0).all())
        np.testing.assert_array_equal(resampled.values, np.round(resampled.values))
        self.assertEqual(resampled.iloc[:, 1].sum(), 0)
        self.assertFalse(resampled.equals(self.table))


if __name__ == '__main__':
    unittest.main()
