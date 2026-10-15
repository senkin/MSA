"""Merging of plot PDFs (common_methods.merge_plots). Run from the repository root:
    python -m unittest discover -s tests -v
"""
import importlib.util
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

os.environ.setdefault('MPLBACKEND', 'Agg')
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
import matplotlib.pyplot as plt  # noqa: E402
from common_methods import merge_plots  # noqa: E402


class MergePlotsTest(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='msa_plots_test_'))
        (self.tmp / 'plots').mkdir()
        for page in (10, 2, 1):
            fig = plt.figure()
            fig.text(0.5, 0.5, f'page {page}')
            fig.savefig(self.tmp / f'plots/sample_{page}.pdf')
            plt.close(fig)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    @unittest.skipUnless(importlib.util.find_spec('pypdf'), 'pypdf is not installed')
    def test_natural_order(self):
        from pypdf import PdfReader
        merged = self.tmp / 'all.pdf'
        merge_plots(str(self.tmp / 'plots/*pdf'), str(merged))
        titles = [page.extract_text().strip() for page in PdfReader(merged).pages]
        self.assertEqual(titles, ['page 1', 'page 2', 'page 10'])

    def test_without_pypdf_warns(self):
        with mock.patch.dict(sys.modules, {'pypdf': None}):
            with self.assertWarnsRegex(UserWarning, 'pypdf is not installed'):
                merge_plots(str(self.tmp / 'plots/*pdf'), str(self.tmp / 'all.pdf'))
        self.assertFalse((self.tmp / 'all.pdf').exists())


if __name__ == '__main__':
    unittest.main()
