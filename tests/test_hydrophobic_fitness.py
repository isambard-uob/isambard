"""Tests the hydrophobic fitness code."""

import unittest
import pathlib

import ampal

from isambard.evaluation import calculate_hydrophobic_fitness


TEST_FILES_PATH = pathlib.Path(__file__).parent / 'testing_files'


class HelixTestCase(unittest.TestCase):
    """Tests for isambard.evaluation.hydrophobic_fitness"""

    def setUp(self):
        self.ctf = ampal.load_pdb(str(TEST_FILES_PATH / '1ctf.pdb'))
        self.ubq = ampal.load_pdb(str(TEST_FILES_PATH / '1ubq.pdb'))
        self.r69 = ampal.load_pdb(str(TEST_FILES_PATH / '1r69.pdb'))
        self.icb = ampal.load_pdb(str(TEST_FILES_PATH / '4icb.pdb'))

    def test_calculate_hf(self):
        """Checks the relative value of HF for 4 test structures."""
        ctf_score = calculate_hydrophobic_fitness(self.ctf)
        ubq_score = calculate_hydrophobic_fitness(self.ubq)
        r69_score = calculate_hydrophobic_fitness(self.r69)
        icb_score = calculate_hydrophobic_fitness(self.icb)
        self.assertTrue(ubq_score < icb_score < r69_score < ctf_score)


__author__ = 'Christopher W. Wood'
