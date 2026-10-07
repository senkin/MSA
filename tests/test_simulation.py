"""Simulated data (bin/simulate_data.py) from attributions that include samples without
mutations, as sparse mutation types (e.g. DBS) often have. Run from the repository root:
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


class EmptySamplesTest(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='msa_simulation_test_'))
        rng = np.random.default_rng(0)
        # unoptimised attribution (mutation counts) of 10 samples, 2 without DBS mutations
        self.activities = pd.DataFrame(rng.integers(5, 50, size=(10, 5)),
                                       columns=[f'DBS{i}' for i in range(1, 6)])
        self.activities.iloc[[3, 7]] = 0

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def simulate(self, activities, seed):
        activities.to_csv(self.tmp / 'activities.csv')
        out = self.tmp / f'out_{seed}'
        env = dict(os.environ, MPLCONFIGDIR=str(self.tmp / '.matplotlib'), MPLBACKEND='Agg')
        args = ['-d', 'SIMX', '-t', 'DBS', '-n', 20, '-p', 'sigTest', '-s', REPO / 'signature_tables',
                '-B', '-z', '--noise_type', 'gaussian', '-Z', 10, '-i', self.tmp / 'activities.csv',
                '-o', out, '--seed', seed]
        result = subprocess.run([sys.executable, str(REPO / 'bin/simulate_data.py'), *map(str, args)],
                                cwd=self.tmp, env=env, capture_output=True, text=True)
        if result.returncode != 0:
            raise AssertionError(f'simulate_data.py failed:\n{result.stdout}\n{result.stderr}')
        counts = pd.read_csv(out / 'SIMX/WGS_SIMX.dinucs.csv', index_col=0)
        weights = pd.read_csv(out / 'SIMX/WGS_SIMX.dinucs.weights.csv', index_col=0)
        return counts, weights, result.stderr

    def simulate_seeds(self, activities):
        # seeds 1 and 8 drew an empty sample last, which used to empty the whole simulated set
        return [self.simulate(activities, seed) for seed in (1, 2, 8)]

    def test_empty_samples_skipped_with_a_warning(self):
        for counts, weights, stderr in self.simulate_seeds(self.activities):
            self.assertIn('2 of 10 input samples have no DBS mutations', stderr)
            self.assertTrue((counts.sum() > 0).all())
            sums = weights.sum(axis=1).values
            self.assertTrue(np.all((sums > 0) & (sums <= 1 + 1e-9)))

    def test_no_empty_samples_no_warning(self):
        counts, _, stderr = self.simulate(self.activities.loc[self.activities.sum(axis=1) > 0], 1)
        self.assertNotIn('no DBS mutations', stderr)
        self.assertTrue((counts.sum() > 0).all())

    def test_all_samples_empty(self):
        # nothing to simulate from: warn and write empty samples instead of failing
        counts, weights, stderr = self.simulate(self.activities * 0, 1)
        self.assertIn('None of the 10 input samples has DBS mutations', stderr)
        self.assertEqual(counts.shape[1], 20)
        self.assertEqual(counts.values.sum(), 0)
        self.assertEqual(weights.values.sum(), 0)


if __name__ == '__main__':
    unittest.main()
