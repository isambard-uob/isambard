import os
import random
import unittest
import warnings

import isambard_dev as isambard
from hypothesis import given, settings
from hypothesis.strategies import integers

warnings.filterwarnings("ignore")


class ForceFieldTestCase(unittest.TestCase):
    """Tests for isambard.buff.BuffForceField."""

    def setUp(self):
        self.ff = isambard.buff.BuffForceField()
        self.topo = isambard.specifications.CoiledCoil(3)
        pdb_path = os.path.join('unit_tests', 'testing_files', '3qy1.pdb')
        self.pdb = isambard.ampal.convert_pdb_to_ampal(pdb_path)

    def test_ff_json(self):
        """Tests if all force fields present are correctly formatted."""
        for ff_type in isambard.buff.force_fields:
            isambard.buff.BuffForceField(force_field=ff_type)

    def test_default_ff(self):
        """Tests if all the parameters are correctly loaded."""
        for res, atoms in self.ff.items():
            if res != "KEY":
                for atom, params in atoms.items():
                    self.assertEqual(len(params), 8)
                    self.assertTrue(type(atom) is str)
                    self.assertTrue(type(params[0]) is str)
                    self.assertTrue(all([(type(x) is int) or (type(x) is float) for x in params[1:]]))

    def test_parameterisation_topo(self):
        """Checks that all atoms that can be parameterised are parameterised."""
        self.topo.assign_force_field(self.ff)
        for atom in self.topo.get_atoms(inc_alt_states=True):
            if atom.element != 'H':
                self.assertTrue(atom._ff_id is not None)

    def test_parameterisation_pdb(self):
        """Checks that all atoms that can be parameterised are parameterised."""
        self.pdb.assign_force_field(self.ff, mol2=True)
        for atom in self.pdb.get_atoms(inc_alt_states=True):
            if atom.element != 'H':
                if atom.ampal_parent.mol_code != 'HOH':
                    self.assertTrue(atom._ff_id is not None)


class InteractionsTestCase(unittest.TestCase):
    """Tests for isambard.buff.find_buff_interactions."""

    def setUp(self):
        self.ff = isambard.buff.BuffForceField()  # Auto loads default
        self.topo = isambard.specifications.CoiledCoil(3)
        self.topo.assign_force_field(self.ff)
        pdb_path = os.path.join('unit_tests', 'testing_files', '3qy1.pdb')
        self.pdb = isambard.ampal.convert_pdb_to_ampal(pdb_path)
        self.pdb.assign_force_field(self.ff)

    def test_basic_intra_format(self):
        """Tests that the interaction tuples are correctly formatted."""
        buff_interactions = isambard.buff.find_intra_ampal(self.topo, self.ff.distance_cutoff)
        for _ in range(100):
            a, b = random.choice(buff_interactions)
            self.assertTrue(type(a) is isambard.ampal.Atom)
            self.assertTrue(type(b) is isambard.ampal.Atom)
            self.assertTrue(a != b)

    def test_basic_inter_format(self):
        """Tests that the interaction tuples are correctly formatted."""
        buff_interactions = isambard.buff.find_inter_ampal(self.topo, self.ff.distance_cutoff)
        for _ in range(100):
            a, b = random.choice(buff_interactions)
            self.assertTrue(type(a) is isambard.ampal.Atom)
            self.assertTrue(type(b) is isambard.ampal.Atom)
            self.assertTrue(a != b)

    @given(integers(min_value=1, max_value=50))
    @settings(max_examples=20)
    def test_intra_num_helix(self, n):
        """Tests that the number of interactions found in a Helix backbone."""
        ia_scaling = lambda x: (((x - 2) * (x - 1)) / 2) * 16 if x > 0.0 else 0.0
        helix = isambard.specifications.Helix(n)
        helix.assign_force_field(self.ff)
        buff_interactions = isambard.buff.find_intra_ampal(helix, 1.52*(n+1))
        self.assertEqual(len(buff_interactions), ia_scaling(n))

    @given(integers(min_value=1, max_value=30), integers(min_value=2, max_value=5))
    @settings(max_examples=10)
    def test_inter_num_cc(self, n, hels):
        """Tests that the number of interactions found between helices in a coiled coil."""
        ia_scaling = lambda x, y: (((x ** 2) * (y - 1) * y) / 2) * 16 if x > 0.0 else 0.0
        cc = isambard.specifications.CoiledCoil.from_parameters(hels, n, 7.0, 180.0, 18.0)
        cc.assign_force_field(self.ff)
        buff_interactions = isambard.buff.find_inter_ampal(cc, 1000)
        self.assertEqual(len(buff_interactions), ia_scaling(n, hels))

    def test_interaction_energy(self):
        """Tests that the interaction energy of a reference structure."""
        buff_score = self.pdb.get_interaction_energy(ff=self.ff)
        self.assertAlmostEqual(buff_score.total_energy, -1005.41, places=2)  # Original scoring function

    def test_internal_energy(self):
        """Tests that the internal energy of a reference structure."""
        buff_score = self.pdb[0].get_internal_energy(ff=self.ff)
        self.assertAlmostEqual(buff_score.total_energy, -3722.49, places=2)  # Original scoring function


if __name__ == '__main__':
    unittest.main()
