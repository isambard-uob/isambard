import os
import random
import unittest

import isambard_dev as isambard


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

    def test_basic_interaction(self):
        """Tests that the interaction tuples are correctly formatted."""
        buff_interactions = isambard.buff.calculate_energy.find_buff_interactions(self.topo, self.ff)
        for _ in range(100):
            a, b = random.choice(buff_interactions)
            self.assertTrue(type(a) is isambard.ampal.Atom)
            self.assertTrue(type(b) is isambard.ampal.Atom)
            self.assertTrue(a != b)

    def test_inter_count_pdb(self):
        """Tests that the number of interactions found in a pdb file."""
        buff_interactions = isambard.buff.calculate_energy.find_buff_interactions(self.pdb, self.ff)
        self.assertEqual(len(buff_interactions), 66546)  # Original scoring function


if __name__ == '__main__':
    unittest.main()
