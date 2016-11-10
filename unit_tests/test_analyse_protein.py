import unittest
import os

import isambard_dev as isambard

_test_file = os.path.join('unit_tests', 'testing_files', '2ebo_1.mmol')
_test_polypeptide = isambard.ampal.convert_pdb_to_ampal(_test_file)[0]


class ResiduesPerTurnTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rpts = list(isambard.analyse_protein.residues_per_turn(p=_test_polypeptide))

    def test_residues_per_turn_length(self):
        self.assertEqual(len(self.rpts), len(_test_polypeptide.primitive))

    def test_residues_per_turn_final_value(self):
        self.assertIsNone(self.rpts[-1])

    def test_residues_per_turn_none_values(self):
        none_values = [x for x in self.rpts if x is None]
        self.assertAlmostEqual(len(none_values), 1)

    def test_residues_per_turn_index_eight(self):
        self.assertAlmostEqual(self.rpts[8], 3.6400862)


class RadiiOfCurvatureTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rocs = list(_test_polypeptide.radii_of_curvature())

    def test_radii_of_curvature_length(self):
        self.assertEqual(len(self.rocs), len(_test_polypeptide.primitive))

    def test_radii_of_curvature_first_and_final_values(self):
        self.assertIsNone(self.rocs[0])
        self.assertIsNone(self.rocs[-1])

    def test_radii_of_curvature_none_values(self):
        none_values = [x for x in self.rocs if x is None]
        self.assertEqual(len(none_values), 2)

    def test_radii_of_curvature_index_eighteen(self):
        self.assertAlmostEqual(self.rocs[18], 17.1549833)


class RisePerResidueTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rprs = list(_test_polypeptide.rise_per_residue())

    def test_rise_per_residue_length(self):
        self.assertEqual(len(self.rprs), len(_test_polypeptide.primitive))

    def test_rise_per_residue_final_value(self):
        self.assertIsNone(self.rprs[-1])

    def test_rise_per_residue_none_values(self):
        none_values = [x for x in self.rprs if x is None]
        self.assertEqual(len(none_values), 1)

    def test_rise_per_residue_index_thirteen(self):
        self.assertAlmostEqual(self.rprs[13], 1.4341967)


__author__ = 'Jack W. Heal'
