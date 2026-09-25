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

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / 'bin'
SIGNATURES = REPO / 'signature_tables'


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


if __name__ == '__main__':
    unittest.main()
