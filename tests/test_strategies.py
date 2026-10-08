"""The signature optimisation strategies of bin/run_NNLS.py: 'addition' and 'add-remove' start
from SBS1 and SBS5 (or the first signature of the table) and add the strong signatures.
Run from the repository root:
    python -m unittest discover -s tests -v
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
import run_NNLS as R  # noqa: E402


def signatures(mutation_type):
    index_col = [0, 1] if mutation_type == 'SBS' else 0
    return pd.read_csv(REPO / f'signature_tables/sigProfiler_{mutation_type}_signatures.csv', index_col=index_col)


def mixture(sigs, weights, burden=5000):
    """Mutation counts of an exact mixture of signatures"""
    weights = pd.Series(weights)
    return (sigs[weights.index] @ weights * burden).values


def optimise(observed, sigs, strategy, weak=0.01, strong=0.01):
    with np.errstate(invalid='ignore'):  # exact fits: Jensen-Shannon takes sqrt of a rounding -0
        return list(R.optimise_signatures(observed, sigs, sigs, strategy=strategy,
                                          weak_threshold=weak, strong_threshold=strong).columns)


class StartSignaturesTest(unittest.TestCase):

    def test_start_signatures(self):
        self.assertEqual(R.addition_start_signatures(['SBS1', 'SBS2', 'SBS5']), ['SBS1', 'SBS5'])
        self.assertEqual(R.addition_start_signatures(['SBS2', 'SBS5', 'SBS13']), ['SBS5'])
        # neither available: the first signature listed
        self.assertEqual(R.addition_start_signatures(['DBS1', 'DBS2']), ['DBS1'])
        self.assertEqual(R.addition_start_signatures([]), [])


class StrategiesTest(unittest.TestCase):

    def test_addition_adds_strong_signatures(self):
        sigs = signatures('SBS')
        observed = mixture(sigs, {'SBS1': 0.2, 'SBS5': 0.3, 'SBS13': 0.5})
        # (started from all 65 signatures before, so nothing was added or removed)
        self.assertEqual(sorted(optimise(observed, sigs, 'addition')), ['SBS1', 'SBS13', 'SBS5'])

    def test_add_remove_drops_weak_start_signatures(self):
        sigs = signatures('SBS')
        observed = mixture(sigs, {'SBS2': 0.6, 'SBS13': 0.4})
        # addition keeps the signatures it starts from, add-remove removes them if weak
        self.assertEqual(sorted(optimise(observed, sigs, 'addition')), ['SBS1', 'SBS13', 'SBS2', 'SBS5'])
        self.assertEqual(sorted(optimise(observed, sigs, 'add-remove')), ['SBS13', 'SBS2'])

    def test_first_signature_without_SBS1_and_SBS5(self):
        sigs = signatures('DBS')
        observed = mixture(sigs, {'DBS2': 0.7, 'DBS4': 0.3}, burden=500)
        self.assertEqual(sorted(optimise(observed, sigs, 'addition')), ['DBS1', 'DBS2', 'DBS4'])
        self.assertEqual(sorted(optimise(observed, sigs, 'add-remove')), ['DBS2', 'DBS4'])

    def test_strong_threshold_stops_the_addition(self):
        sigs = signatures('SBS')
        observed = mixture(sigs, {'SBS1': 0.45, 'SBS5': 0.5, 'SBS13': 0.05})
        # SBS13 improves the similarity by less than the strong threshold
        self.assertEqual(sorted(optimise(observed, sigs, 'addition', strong=0.5)), ['SBS1', 'SBS5'])
        self.assertIn('SBS13', optimise(observed, sigs, 'addition', strong=0.001))


if __name__ == '__main__':
    unittest.main()
