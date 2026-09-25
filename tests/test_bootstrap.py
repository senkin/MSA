"""Bootstrap resampling (run_NNLS.py --bootstrap_method) on every mutation type and SBS context.

Run from the repository root:
    python -m unittest discover -s tests -v
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / 'bin'
INPUTS = REPO / 'input_mutation_tables/SIM_test'
SIGNATURES = REPO / 'signature_tables'
sys.path.insert(0, str(BIN))
import run_NNLS as R  # noqa: E402

# (mutation type, context, input table suffix, index columns, signature prefix)
DATASETS = [
    ('SBS', 96, '96', [0, 1], 'sigTest'),
    ('SBS', 192, '192', [0, 1, 2], 'sigTest'),
    ('SBS', 288, '288', [0, 1, 2], 'sigTest'),
    ('SBS', 1536, '1536', 0, 'sigTest'),
    ('SBS', 4608, '4608', 0, 'sigTest'),
    ('DBS', 96, 'dinucs', 0, 'sigTest'),
    ('ID', 96, 'indels', 0, 'sigTest'),
    ('SV', 96, 'SV', 0, 'sigProfiler'),
    ('CNV', 96, 'CNV', 0, 'sigProfiler'),
]
WEAK, STRONG = [0.0, 0.01], [0.0]


def label(dataset):
    mutation_type, context = dataset[:2]
    return f'SBS{context}' if mutation_type == 'SBS' else mutation_type


def load(mutation_type, context, suffix, index_col, prefix):
    """(mutations, signatures) of SIM_test, with a sample without mutations in column 1."""
    name = f'{prefix}_{mutation_type}' + (f'_{context}' if mutation_type == 'SBS' and context != 96 else '')
    signatures = pd.read_csv(SIGNATURES / f'{name}_signatures.csv', index_col=index_col)
    mutations = pd.read_csv(INPUTS / f'WGS_SIM_test.{suffix}.csv', index_col=index_col)
    mutations.iloc[:, 1] = 0
    return mutations, signatures


def ignore_empty_sample_warnings(test):
    """For the rest of the test, silence the warning run_NNLS.py gives each time it skips a
    sample without mutations (column 1 here, or a resample that drew no mutations)."""
    filters = warnings.catch_warnings()
    filters.__enter__()
    test.addCleanup(filters.__exit__, None, None, None)
    warnings.filterwarnings('ignore', 'Sample .*: Zero mutations', UserWarning)


def make_args(method):
    return argparse.Namespace(optimise_signatures=True, optimisation_strategy='removal',
                              weak_threshold=WEAK[0], strong_threshold=STRONG[0],
                              normalise_mutations=False, verbose=False, use_gpu=False,
                              gpu_batch_size=0, threshold_scan='trajectory', bootstrap_method=method)


def central_fit(mutations, signatures, method, weak=WEAK, strong=STRONG):
    """The central attributions bootstrap_residuals resamples (None for the other methods)."""
    if method != 'bootstrap_residuals':
        return None
    return R.attribute_at_thresholds(mutations.reindex(signatures.index), signatures,
                                     list(range(signatures.shape[1])), make_args(method), weak, strong)[0]


def resample(mutations, signatures, method, seed, rows=None):
    np.random.seed(seed)
    if method != 'bootstrap_residuals':
        return R.bootstrap_mutation_table(mutations, method)
    _, _, _, residuals, fitted = central_fit(mutations, signatures, method)[(WEAK[0], STRONG[0])]
    return R.bootstrap_mutation_table(None, method, fitted=fitted, residuals=residuals, rows=rows)


class ResamplingTest(unittest.TestCase):
    """bootstrap_mutation_table(): one resampled table per method."""

    @classmethod
    def setUpClass(cls):
        cls.data = {label(d): load(*d) for d in DATASETS}

    def setUp(self):
        ignore_empty_sample_warnings(self)

    def test_every_method(self):
        for name, (mutations, signatures) in self.data.items():
            for method in R.BOOTSTRAP_METHODS:
                with self.subTest(dataset=name, method=method):
                    resampled = resample(mutations, signatures, method, seed=0)
                    self.assertTrue(resampled.index.equals(mutations.index))
                    self.assertTrue(resampled.index.is_unique)
                    self.assertTrue(resampled.columns.equals(mutations.columns))
                    values = resampled.values.astype(float)
                    self.assertTrue(np.isfinite(values).all() and (values >= 0).all())
                    self.assertEqual(values[:, 1].sum(), 0)  # nothing to resample
                    self.assertFalse(np.array_equal(values, mutations.values))
                    if method != 'bootstrap_residuals':
                        np.testing.assert_array_equal(values, np.round(values))
                    # reproducible from the seed alone, and a new seed gives a new resample
                    np.testing.assert_array_equal(resample(mutations, signatures, method, seed=0).values,
                                                  resampled.values)
                    self.assertFalse(np.array_equal(resample(mutations, signatures, method, seed=1).values,
                                                    resampled.values))

    def test_multinomial_keeps_totals(self):
        for name, (mutations, signatures) in self.data.items():
            with self.subTest(dataset=name):
                resampled = resample(mutations, signatures, 'multinomial', seed=0)
                np.testing.assert_array_equal(resampled.sum().values, mutations.sum().values)

    def test_poisson_keeps_empty_channels(self):
        for name, (mutations, signatures) in self.data.items():
            with self.subTest(dataset=name):
                resampled = resample(mutations, signatures, 'poisson', seed=0)
                self.assertTrue((resampled.values[mutations.values == 0] == 0).all())

    def test_classic_resamples_channels(self):
        # whole channels (rows) drawn with replacement; the channel labels stay in place
        for name, (mutations, signatures) in self.data.items():
            with self.subTest(dataset=name):
                resampled = resample(mutations, signatures, 'classic', seed=0)
                self.assertLessEqual(set(map(tuple, resampled.values)), set(map(tuple, mutations.values)))
                rows = np.arange(len(mutations))[::-1]
                np.testing.assert_array_equal(R.bootstrap_mutation_table(mutations, 'classic', rows=rows).values,
                                              mutations.values[rows])

    def test_multinomial_weight_resamples_within_samples(self):
        for name, (mutations, signatures) in self.data.items():
            with self.subTest(dataset=name):
                resampled = resample(mutations, signatures, 'multinomial_weight', seed=0)
                for sample in mutations.columns:
                    self.assertLessEqual(set(resampled[sample]), set(mutations[sample]))

    def test_residuals_added_to_the_central_fit(self):
        for name, (mutations, signatures) in self.data.items():
            with self.subTest(dataset=name):
                _, _, _, residuals, fitted = central_fit(mutations, signatures,
                                                         'bootstrap_residuals')[(WEAK[0], STRONG[0])]
                observed = mutations.reindex(signatures.index)
                # without resampling, the fit plus its residuals is the observed table
                identity = np.arange(len(observed))
                unchanged = R.bootstrap_mutation_table(None, 'bootstrap_residuals', fitted=fitted,
                                                       residuals=residuals, rows=identity)
                np.testing.assert_allclose(unchanged.values, observed.values, atol=1e-9)
                self.assertTrue(unchanged.index.equals(observed.index))
                rows = np.random.default_rng(0).integers(len(observed), size=len(observed))
                resampled = R.bootstrap_mutation_table(None, 'bootstrap_residuals', fitted=fitted,
                                                       residuals=residuals, rows=rows)
                np.testing.assert_allclose(resampled.values,
                                           np.clip(fitted.values + residuals.values[rows], 0, None))

    def test_invalid_requests_raise(self):
        mutations, _ = self.data['SBS96']
        with self.assertRaises(ValueError):
            R.bootstrap_mutation_table(mutations, 'bootstrap_residuals')
        with self.assertRaises(ValueError):
            R.bootstrap_mutation_table(mutations, 'unknown')


class BootstrapAttributionTest(unittest.TestCase):
    """bootstrap_attributions(): the attributions of one resample, at several thresholds."""

    @classmethod
    def setUpClass(cls):
        cls.data = {label(d): load(*d) for d in DATASETS}

    def setUp(self):
        ignore_empty_sample_warnings(self)

    def attribute(self, mutations, signatures, method, weak, seed):
        central = central_fit(mutations, signatures, method, weak)
        np.random.seed(seed)
        return R.bootstrap_attributions(mutations, signatures, list(range(signatures.shape[1])),
                                        make_args(method), weak, STRONG, central)

    def test_every_method(self):
        for name, (mutations, signatures) in self.data.items():
            for method in R.BOOTSTRAP_METHODS:
                with self.subTest(dataset=name, method=method):
                    results = self.attribute(mutations, signatures, method, WEAK, seed=3)
                    self.assertEqual(set(results), {(w, s) for w in WEAK for s in STRONG})
                    for (weak, strong), outputs in results.items():
                        weights = outputs[0]
                        self.assertEqual(list(weights.index), list(mutations.columns))
                        self.assertEqual(list(weights.columns), list(signatures.columns))
                        self.assertTrue(np.isfinite(weights.values).all() and (weights.values >= 0).all())
                        sums = weights.sum(axis=1).values
                        np.testing.assert_allclose(sums[sums > 0], 1)
                        self.assertEqual(sums[1], 0)
                        # every threshold sees the same resample as a run at that threshold alone
                        alone = self.attribute(mutations, signatures, method, [weak], seed=3)[(weak, strong)]
                        for got, expected in zip(outputs[:3], alone[:3]):
                            np.testing.assert_allclose(got.values, expected.values, rtol=0, atol=1e-12)

    def test_empty_sample_is_skipped_with_a_warning(self):
        mutations, signatures = self.data['SBS96']
        with self.assertWarnsRegex(UserWarning, 'Sample 1: Zero mutations'):
            self.attribute(mutations, signatures, R.BOOTSTRAP_METHODS[0], WEAK, seed=0)

    def test_unresampled_residuals_give_the_central_attribution(self):
        method = 'bootstrap_residuals'
        for name, (mutations, signatures) in self.data.items():
            with self.subTest(dataset=name):
                central = central_fit(mutations, signatures, method)
                # each threshold pair adds back the residuals of its own central fit
                with mock.patch.object(R, 'resample_channels', np.arange):
                    results = R.bootstrap_attributions(mutations, signatures, list(range(signatures.shape[1])),
                                                       make_args(method), WEAK, STRONG, central)
                for pair, outputs in results.items():
                    np.testing.assert_allclose(outputs[0].values, central[pair][0].values, atol=1e-9)


def run_nnls(*args, cwd):
    env = dict(os.environ, MPLCONFIGDIR=str(Path(cwd) / '.matplotlib'), MPLBACKEND='Agg')
    result = subprocess.run([sys.executable, str(BIN / 'run_NNLS.py'), *map(str, args)],
                            cwd=cwd, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError(f'run_NNLS.py failed:\n{result.stdout}\n{result.stderr}')


class CommandLineTest(unittest.TestCase):
    """run_NNLS.py as the bootstrap processes run it (-B --n_bootstrap, one folder per threshold)."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='msa_bootstrap_test_'))
        shutil.copytree(INPUTS, self.tmp / 'inputs/SIM_test')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def bootstrap(self, dataset, out, *options):
        mutation_type, context, _, _, prefix = dataset
        run_nnls('-B', '-d', 'SIM_test', '-t', mutation_type, '-c', context, '-x', '-p', prefix,
                 '-i', self.tmp / 'inputs', '-s', SIGNATURES, '-o', out, '--seed', 1, *options, cwd=self.tmp)

    def check_outputs(self, dataset, out, n_bootstrap):
        mutation_type, _, _, _, prefix = dataset
        mutations, signatures = load(*dataset)
        for weak in WEAK:
            for strong in STRONG:
                pair = f'{R.format_threshold(weak)}_{R.format_threshold(strong)}'
                folder = Path(out) / f'SIM_test_{dataset[1]}_NNLS_{pair}/bootstrap_output'
                self.assertEqual(len(list(folder.iterdir())), 3 * n_bootstrap)
                for i in range(1, n_bootstrap + 1):
                    weights = pd.read_csv(folder / f'output_SIM_test_{mutation_type}_{pair}_{i}_weights_table.csv',
                                          index_col=0)
                    self.assertEqual(list(weights.index.astype(str)), list(mutations.columns))
                    self.assertEqual(list(weights.columns), list(signatures.columns))
                    self.assertFalse(weights.isna().any().any())

    def test_default_bootstrap_every_mutation_type(self):
        for dataset in DATASETS:
            with self.subTest(dataset=label(dataset)):
                out = self.tmp / label(dataset)
                self.bootstrap(dataset, out, '--n_bootstrap', 2, '-W', *WEAK, '-S', *STRONG, '--add_suffix')
                self.check_outputs(dataset, out, 2)

    def test_every_method(self):
        for method in R.BOOTSTRAP_METHODS:
            if method == 'binomial':
                continue  # the default, run above
            with self.subTest(method=method):
                out = self.tmp / method
                self.bootstrap(DATASETS[0], out, '--n_bootstrap', 2, '--bootstrap_method', method,
                               '-W', *WEAK, '-S', *STRONG, '--add_suffix')
                self.check_outputs(DATASETS[0], out, 2)

    def test_single_bootstrap_is_an_iteration_of_the_loop(self):
        # -B alone runs iteration --bootstrap_start_index of the in-process loop, also with -n
        for method in ('binomial', 'bootstrap_residuals'):
            with self.subTest(method=method):
                loop, single = self.tmp / f'loop_{method}', self.tmp / f'single_{method}'
                options = ('--bootstrap_method', method, '-W', 0.01, '-S', 0, '-n', 3)
                self.bootstrap(DATASETS[1], loop, '--n_bootstrap', 2, '--bootstrap_start_index', 4, *options)
                self.bootstrap(DATASETS[1], single, '--bootstrap_start_index', 5, *options)
                for table in ('weights_table', 'stat_info'):
                    self.assertEqual(
                        (loop / f'SIM_test/bootstrap_output/output_SIM_test_SBS_5_{table}.csv').read_bytes(),
                        (single / f'SIM_test/output_SIM_test_SBS_{table}.csv').read_bytes())

    def test_central_run_leaves_inputs_untouched(self):
        before = sorted(os.listdir(self.tmp / 'inputs/SIM_test'))
        run_nnls('-d', 'SIM_test', '-t', 'SBS', '-c', 96, '-x', '-W', 0.01, '-S', 0, '-p', 'sigTest',
                 '-i', self.tmp / 'inputs', '-s', SIGNATURES, '-o', self.tmp / 'central', cwd=self.tmp)
        self.assertEqual(sorted(os.listdir(self.tmp / 'inputs/SIM_test')), before)
        self.assertTrue((self.tmp / 'central/SIM_test/output_SIM_test_SBS_residuals.csv').exists())


if __name__ == '__main__':
    unittest.main()
