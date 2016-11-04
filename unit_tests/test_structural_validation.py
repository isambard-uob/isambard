import unittest
import itertools
import os
import random
import numpy

import isambard_dev as isambard
from unit_tests.random_isambard_objects import random_angles, random_floats


class PolypeptideStructuralValidationTestCase(unittest.TestCase):
    """ Tests for the valid_backbone_distance method of Polypeptide class. """
    def setUp(self):
        test_files = [os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', x)
                      for x in ['1ek9.pdb', '2ht0.pdb', '3qy1.pdb']]
        test_structures = [isambard.ampal.convert_pdb_to_ampal(x) for x in test_files]
        self.test_polypeptides = [
            p for p in itertools.chain(*test_structures) if type(p) == isambard.ampal.Polypeptide]

    def test_valid_distances_testing_files(self):
        # checks for valid bond distances in testing files.
        valid_distances = [p.valid_backbone_bond_lengths() for p in self.test_polypeptides]
        self.assertTrue(all(valid_distances))

    def test_invalid_distance_natural(self):
        polypeptides = self.test_polypeptides[:]
        # translate random monomer
        for p in polypeptides:
            p[random.choice(range(1, len(p)))].translate([0, 0, 10])
            self.assertFalse(p.valid_backbone_bond_lengths())

    def test_valid_angles_testing_files(self):
        # checks for valid bond angles in testing files.
        valid_angles = [p.valid_backbone_bond_angles() for p in self.test_polypeptides]
        self.assertTrue(all(valid_angles))

    def test_invalid_helix(self):
        h1 = isambard.specifications.Helix()
        h2 = isambard.specifications.Helix()
        h1.extend(h2)
        self.assertFalse(h1.valid_backbone_bond_lengths())
        self.assertFalse(h1.valid_backbone_bond_angles())


class PolypeptideJoinTestCase(unittest.TestCase):
    def setUp(self):
        self.num_tests = 10
        self.torsion_angles = [random_angles(n=3, min_val=-179, max_val=180) for _ in range(self.num_tests)]
        self.peptide_bond_lengths = random_floats(n=self.num_tests, min_val=1, max_val=3)

    def test_c_join_helices(self):
        tests = []
        for omega, phi, psi in self.torsion_angles:
            len_h1 = 3
            h1 = isambard.specifications.Helix(aa=len_h1)
            h2 = isambard.specifications.Helix(aa=len_h1)
            h1.c_join(h2, psi=psi, phi=phi, omega=omega)
            measured_torsion_angles = isambard.analyse_protein.measure_torsion_angles(h1)
            measured_omega, measured_phi = measured_torsion_angles[len_h1][:2]
            measured_psi = measured_torsion_angles[len_h1 - 1][2]
            tests.append(numpy.allclose([measured_omega, measured_phi, measured_psi], [omega, phi, psi]))
        self.assertTrue(all(tests))

    def test_n_join_helices(self):
        tests = []
        for omega, phi, psi in self.torsion_angles:
            len_h1 = 3
            h1 = isambard.specifications.Helix(aa=len_h1)
            h2 = isambard.specifications.Helix(aa=len_h1)
            h1.n_join(h2, psi=psi, phi=phi, omega=omega)
            measured_torsion_angles = isambard.analyse_protein.measure_torsion_angles(h1)
            measured_omega, measured_phi = measured_torsion_angles[len_h1][:2]
            measured_psi = measured_torsion_angles[len_h1 - 1][2]
            tests.append(numpy.allclose([measured_omega, measured_phi, measured_psi], [omega, phi, psi]))
        self.assertTrue(all(tests))

    def test_c_join_residue_recursive(self):
        tests = []
        h = isambard.specifications.Helix(aa=1)
        for omega, phi, psi in self.torsion_angles:
            h2 = isambard.specifications.Helix(aa=1)
            h.c_join(h2, psi=psi, phi=phi, omega=omega)
            measured_torsion_angles = isambard.analyse_protein.measure_torsion_angles(h)
            measured_omega, measured_phi = measured_torsion_angles[-1][:2]
            measured_psi = measured_torsion_angles[-2][2]
            tests.append(numpy.allclose([measured_omega, measured_phi, measured_psi], [omega, phi, psi]))
        self.assertTrue(all(tests))

    def test_n_join_residue_recursive(self):
        tests = []
        h = isambard.specifications.Helix(aa=1)
        for omega, phi, psi in self.torsion_angles:
            h2 = isambard.specifications.Helix(aa=1)
            h.n_join(h2, psi=psi, phi=phi, omega=omega)
            measured_torsion_angles = isambard.analyse_protein.measure_torsion_angles(h)
            measured_omega, measured_phi = measured_torsion_angles[1][:2]
            measured_psi = measured_torsion_angles[0][2]
            tests.append(numpy.allclose([measured_omega, measured_phi, measured_psi], [omega, phi, psi]))
        self.assertTrue(all(tests))

    def test_c_join_bond_length(self):
        bond_lengths = []
        for i in range(self.num_tests):
            omega, phi, psi = self.torsion_angles[i]
            h1 = isambard.specifications.Helix(aa=1)
            h2 = isambard.specifications.Helix(aa=1)
            h1.c_join(h2, psi=psi, phi=phi, omega=omega, c_n_length=self.peptide_bond_lengths[i])
            bond_lengths.append(isambard.geometry.distance(h1[0]['C'], h1[1]['N']))
        at = numpy.allclose(bond_lengths, self.peptide_bond_lengths)
        self.assertTrue(at)

    def test_n_join_bond_length(self):
        bond_lengths = []
        for i in range(self.num_tests):
            omega, phi, psi = self.torsion_angles[i]
            h1 = isambard.specifications.Helix(aa=1)
            h2 = isambard.specifications.Helix(aa=1)
            h1.n_join(h2, psi=psi, phi=phi, omega=omega, c_n_length=self.peptide_bond_lengths[i])
            bond_lengths.append(isambard.geometry.distance(h1[0]['C'], h1[1]['N']))
        at = numpy.allclose(bond_lengths, self.peptide_bond_lengths)
        self.assertTrue(at)

    def test_valid_bond_lengths_and_angles_c_join(self):
        h = isambard.specifications.Helix(aa=1)
        for omega, phi, psi in self.torsion_angles:
            h2 = isambard.specifications.Helix(aa=1)
            h.c_join(h2, psi=psi, phi=phi, omega=omega)
            self.assertTrue(h.valid_backbone_bond_angles())
            self.assertTrue(h.valid_backbone_bond_lengths())

    def test_valid_bond_lengths_and_angles_n_join(self):
        h = isambard.specifications.Helix(aa=1)
        for omega, phi, psi in self.torsion_angles:
            h2 = isambard.specifications.Helix(aa=1)
            h.n_join(h2, psi=psi, phi=phi, omega=omega)
            self.assertTrue(h.valid_backbone_bond_angles())
            self.assertTrue(h.valid_backbone_bond_lengths())

__author__ = 'Jack W. Heal'
