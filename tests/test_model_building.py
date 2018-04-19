import random
import unittest
import numpy
import pathlib

import ampal

import isambard
import isambard.specifications


TEST_FILES_PATH = pathlib.Path(__file__).parent / 'testing_files'


def random_angles(n=1, min_val=0, max_val=180, radians=False):
    angles = [(random.random() * random.choice(range(abs(max_val - min_val)))
               ) + min_val for _ in range(n)]
    if radians:
        angles = [numpy.rad2deg(x) for x in angles]
    return angles


def random_vectors(n=1, min_val=-100, max_val=100, vector_length=3):
    return [[(random.random() * random.choice(range(abs(max_val - min_val)))) + min_val
             for _ in range(vector_length)] for _ in range(n)]


def random_integer_vectors(n=1, min_val=-100, max_val=100, vector_length=3):
    return [[random.choice(range(min_val, max_val)) for _ in range(vector_length)]
            for _ in range(n)]


def random_floats(n=1, min_val=-100, max_val=100):
    return [(random.random() * random.choice(range(abs(max_val - min_val)))) + min_val
            for _ in range(n)]


def random_helical_helices(n=1, min_residues=10, max_residues=20):
    # build random HelicalHelix objects
    major_radii = random_floats(n=n, min_val=1, max_val=100)
    major_pitches = random_floats(n=n, min_val=100, max_val=1000)
    aas = [random.choice(range(min_residues, max_residues)) for _ in range(n)]
    phi_c_alphas = random_angles(n=n, min_val=-179, max_val=180)
    orientations = [random.choice([-1, 1]) for _ in range(n)]
    # minor_repeat can't be set to zero - must be set to None.
    minor_repeats = random_floats(n=n, min_val=0, max_val=100)
    zero_indices = [i for i, x in enumerate(minor_repeats) if x == 0.0]
    for i in zero_indices:
        minor_repeats[i] = None
    minor_helix_types = [random.choice(
        ['alpha', 'pi', 'PPII', 'collagen']) for _ in range(n)]
    major_handedness = [random.choice(['l', 'r']) for _ in range(n)]
    hhs = [isambard.specifications.helix.HelicalHelix(
        aa=aas[i], major_pitch=major_pitches[i],
        major_radius=major_radii[i], major_handedness=major_handedness[i],
        minor_helix_type=minor_helix_types[
            i], orientation=orientations[i],
        phi_c_alpha=phi_c_alphas[i], minor_repeat=minor_repeats[i])
        for i in range(n)]
    return hhs


def check_dihedrals(polypeptide, phi=-65.0, omega=-178.0, psi=-41.0, atol=1):
    """ Returns True if all backbone dihedrals of the polypeptide are closer than atol degrees to the provided values.

    Notes
    -----
    Default values are the ideal dihedral values for the alpha helix.

    """
    omegas, phis, psis = zip(
        *ampal.analyse_protein.measure_torsion_angles(polypeptide))
    omega_check = numpy.allclose(
        omegas[1:], [omega] * (len(polypeptide) - 1), atol=atol)
    phi_check = numpy.allclose(
        phis[1:], [phi] * (len(polypeptide) - 1), atol=atol)
    psi_check = numpy.allclose(
        psis[:-1], [psi] * (len(polypeptide) - 1), atol=atol)
    return all([omega_check, phi_check, psi_check])


class HelixTestCase(unittest.TestCase):
    """Tests for isambard.tools.tool_geometry.dihedral"""

    def setUp(self):
        n = 50
        self.num_tests = n
        self.v1s = random_vectors(n)
        self.v2s = random_vectors(n)
        self.int_v1s = random_integer_vectors(n=n)
        self.int_v2s = random_integer_vectors(n=n)
        self.angles = random_angles(n)
        self.helix_lengths = [random.choice(range(1, 101)) for _ in range(n)]
        self.translations = [random_floats(n)] * 3
        self.rotations = [random_angles(n)] * 3
        self.backbone_angle_tolerance = 5

    def test_alpha_simple(self):
        test_helix = isambard.specifications.Helix(aa=10)
        self.assertEqual(len(test_helix), 10)
        self.assertTrue(test_helix.valid_backbone_bond_lengths())
        self.assertTrue(test_helix.valid_backbone_bond_angles())
        self.assertTrue(check_dihedrals(test_helix))

    def test_alpha_off_z(self):
        test_helix = isambard.specifications.Helix.from_start_and_end(
            (12.3, -1.12, 8.2), (-1.1, 20.2, -12.2), aa=30)
        self.assertEqual(len(test_helix), 30)
        self.assertTrue(test_helix.valid_backbone_bond_lengths())
        self.assertTrue(test_helix.valid_backbone_bond_angles())
        self.assertTrue(check_dihedrals(test_helix))

    def test_random_length_from_o(self):
        """Build a helix of a random length in a random direction from the origin."""
        for i, v in enumerate(self.v1s):
            h_len = self.helix_lengths[i]
            test_helix = isambard.specifications.Helix.from_start_and_end(
                (0, 0, 0), v, aa=h_len)
            reverse_test_helix = isambard.specifications.Helix.from_start_and_end(
                v, (0, 0, 0), aa=h_len)
            self.assertEqual(len(test_helix), h_len)
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))
            self.assertEqual(len(reverse_test_helix), h_len)
            self.assertTrue(reverse_test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(reverse_test_helix))

    def test_random_length_direction(self):
        """Build a helix of a random length in a random position."""
        for i in range(len(self.v1s)):
            h_len = self.helix_lengths[i]
            v1 = self.int_v1s[i]
            v2 = self.int_v2s[i]
            test_helix = isambard.specifications.Helix.from_start_and_end(
                v1, v2, aa=h_len)
            self.assertEqual(len(test_helix), h_len)
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))

    def test_random_length_direction_integer_vectors(self):
        """Build a helix of a random length in a random position."""
        for i in range(len(self.v1s)):
            h_len = self.helix_lengths[i]
            v1 = self.v1s[i]
            v2 = self.v2s[i]
            test_helix = isambard.specifications.Helix.from_start_and_end(
                v1, v2, aa=h_len)
            self.assertEqual(len(test_helix), h_len)
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))

    def test_alpha_trans(self):
        test_helix = isambard.specifications.Helix(aa=30)
        self.assertEqual(len(test_helix), 30)
        for i in range(self.num_tests):
            test_helix.translate(
                test_helix.axis.unit_tangent * self.translations[0][i])
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))
            test_helix.translate(
                test_helix.axis.unit_normal * self.translations[1][i])
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))
            test_helix.translate(
                test_helix.axis.unit_binormal * self.translations[2][i])
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))

    def test_alpha_rot(self):
        test_helix = isambard.specifications.Helix(aa=30)
        for i in range(self.num_tests):
            test_helix.rotate(angle=self.rotations[0][i],
                              axis=test_helix.axis.unit_tangent,
                              point=test_helix.helix_start)
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))
            test_helix.rotate(angle=self.rotations[1][i],
                              axis=test_helix.axis.unit_normal,
                              point=test_helix.helix_start)
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))
            test_helix.rotate(angle=self.rotations[2][i],
                              axis=test_helix.axis.unit_tangent,
                              point=test_helix.helix_start)
            self.assertTrue(test_helix.valid_backbone_bond_lengths())
            self.assertTrue(test_helix.valid_backbone_bond_angles(
                atol=self.backbone_angle_tolerance))
            self.assertTrue(check_dihedrals(test_helix))

    def test_pi_simple(self):
        test_helix = isambard.specifications.Helix(aa=30, helix_type='pi')
        self.assertEqual(len(test_helix), 30)
        self.assertTrue(test_helix.valid_backbone_bond_lengths())
        self.assertTrue(test_helix.valid_backbone_bond_angles(
            atol=self.backbone_angle_tolerance))
        self.assertTrue(check_dihedrals(
            test_helix, phi=-75.7, omega=-171.8, psi=-54.8))


class HelicalHelixTestCase(unittest.TestCase):
    """Tests for HelicalHelix sse."""

    def setUp(self):
        n = 25
        self.helical_helices = random_helical_helices(n=n)

    def test_phi_c_alpha(self):
        """ Test that phi_c_alpha used to build is the same as that calculated using get_orient_angle """
        # TODO: THIS TEST_HELICES CONDITION SHOULD BE REMOVED ONCE THE AP BUILDING HAS BE REVISED
        test_helices = [x for x in self.helical_helices if x.orientation == 1]
        calculated_phi_c_alphas = [h.get_orient_angle() for h in test_helices]
        actual_phi_c_alphas = [h.phi_c_alpha for h in test_helices]
        at = numpy.allclose(calculated_phi_c_alphas, actual_phi_c_alphas)
        self.assertTrue(at)


class CoiledCoilTestCase(unittest.TestCase):
    """Tests for the CoiledCoil specifications."""

    def setUp(self):
        n = 25
        self.helical_helices = random_helical_helices(n=n)

    def test_from_polymers(self):
        for _ in range(10):
            oligomer_state = random.choice(range(1, 5))
            coiled_coil = isambard.specifications.CoiledCoil.from_polymers(
                random.sample(self.helical_helices, oligomer_state))
            tests = []
            ''.join(coiled_coil.major_handedness) == ''.join(
                [x.major_handedness for x in coiled_coil._molecules])
            tests.append(numpy.allclose(coiled_coil.major_radii, [
                         x.major_radius for x in coiled_coil._molecules]))
            tests.append(numpy.allclose(coiled_coil.major_pitches, [
                         x.major_pitch for x in coiled_coil._molecules]))
            self.assertTrue(all(tests))


class TABTestCase(unittest.TestCase):

    def setUp(self):
        self.cis_tas = [[-179, 120, -40], [0, -60, 20]]
        self.cis_dipeptide = isambard.specifications.TAPolypeptide(
            self.cis_tas)
        self.trans_tas = [[0, -60, 20], [-179, 120, -40]]
        self.trans_dipeptide = isambard.specifications.TAPolypeptide(
            self.trans_tas)
        test_file = str(TEST_FILES_PATH / '1ek9.pdb')
        test_structure = ampal.load_pdb(test_file)
        self.test_polypeptides = [
            p for p in test_structure if isinstance(p, ampal.Polypeptide)]

    def test_tapolypeptide(self):
        """Testing the build accuracy for TAPolypeptide."""
        torsion_angles = [(-178, -65.0, -41.0)] * 10
        test_pp = isambard.specifications.TAPolypeptide(torsion_angles)
        self.assertEqual(len(torsion_angles), len(test_pp))
        self.assertTrue(test_pp.valid_backbone_bond_lengths(atol=0.01))
        self.assertTrue(test_pp.valid_backbone_bond_angles(atol=5))

    def test_random_ta_polypeptide(self):
        """Testing the build accuracy for a random TAPolypeptide."""
        torsion_angles = [
            (random.choice(range(-180, 180)),
             random.choice(range(-180, 180)),
             random.choice(range(-180, 180))) for _ in range(100)]
        test_pp = isambard.specifications.TAPolypeptide(torsion_angles)
        self.assertEqual(len(torsion_angles), len(test_pp))
        self.assertTrue(test_pp.valid_backbone_bond_lengths(atol=0.01))
        self.assertTrue(test_pp.valid_backbone_bond_angles(atol=5))

    def test_cis_dipeptide_torsion_angles(self):
        self.assertTrue(check_dihedrals(
            self.cis_dipeptide, psi=self.cis_tas[0][2],
            omega=self.cis_tas[1][0], phi=self.cis_tas[1][1], atol=0.0001))

    def test_trans_dipeptide_torsion_angles(self):
        self.assertTrue(check_dihedrals(
            self.trans_dipeptide, psi=self.trans_tas[0][2],
            omega=self.trans_tas[1][0], phi=self.trans_tas[1][1], atol=0.0001))

    def test_cis_didpeptide_valid_backbone(self):
        """ Strict tests for backbone bond angles and distances for cis dipeptide. """
        self.assertTrue(
            self.cis_dipeptide.valid_backbone_bond_angles(atol=0.001))
        self.assertTrue(
            self.cis_dipeptide.valid_backbone_bond_lengths(atol=0.00001))

    def test_trans_didpeptide_valid_backbone(self):
        """ Strict tests for backbone bond angles and distances for cis dipeptide. """
        self.assertTrue(
            self.trans_dipeptide.valid_backbone_bond_angles(atol=0.001))
        self.assertTrue(
            self.trans_dipeptide.valid_backbone_bond_lengths(atol=0.00001))

    def test_from_polypeptide_angles_and_lengths(self):
        for p in self.test_polypeptides:
            tap = isambard.specifications.TAPolypeptide.from_polypeptide(p)
            for k, v in p.backbone_bond_angles.items():
                self.assertTrue(numpy.allclose(tap.backbone_bond_angles[k], v))
            for k, v in p.backbone_bond_lengths.items():
                self.assertTrue(numpy.allclose(
                    tap.backbone_bond_lengths[k], v))

    def test_from_polypeptide_rmsd(self):
        for p in self.test_polypeptides:
            tap = isambard.specifications.TAPolypeptide.from_polypeptide(p)
            p_coords = [x._vector for x in p.backbone.get_atoms()]
            tap_coords = [x._vector for x in tap.get_atoms()]
            self.assertTrue(numpy.isclose(ampal.geometry.rmsd(
                p_coords, tap_coords), 0.0, atol=0.001))


__author__ = 'Christopher W. Wood, Jack W. Heal'
