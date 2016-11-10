from collections import Counter
import os
import unittest

import isambard_dev as isambard


class PdbParserTestCase(unittest.TestCase):
    """Tests for isambard.ampal.PdbParser"""

    def check_ampal_contents(self, pdb_path):
        """Tests if all atoms present are correctly formatted."""
        with open(pdb_path, 'r') as inf:
            pdb_lines = inf.readlines()
        pdb_atom_numbers = set()
        pdb_monomer_labels = {}
        pdb_has_alt_conf = set()
        for line in pdb_lines:
            record_name = line[:6].strip()
            if (record_name == 'ATOM') or (record_name == 'HETATM'):
                atom_number = int(line[6:11].strip())
                pdb_atom_numbers.add(atom_number)

                monomer_number = int(line[22:26].strip())
                chain = line[21].strip()
                monomer_type = line[17:20].strip()
                pdb_monomer_labels[(chain, monomer_number)] = monomer_type

                alt_loc = line[16].strip()
                if alt_loc:
                    pdb_has_alt_conf.add((chain, monomer_number))

        ampal = isambard.ampal.convert_pdb_to_ampal(pdb_path)
        ampal_atom_numbers = set()
        ampal_has_alt_conf = set()
        for atom in ampal.get_atoms(inc_alt_states=True):
            ampal_atom_numbers.add(atom.id)

        ampal_monomer_labels = []
        for monomer in ampal.get_monomers():
            ampal_monomer_labels.append(monomer.mol_code)
            if len(monomer.states) > 1:
                ampal_has_alt_conf.add((monomer.ampal_parent.id, int(monomer.id)))

        self.assertEqual(pdb_atom_numbers, ampal_atom_numbers)
        self.assertEqual(Counter(pdb_monomer_labels.values()), Counter(ampal_monomer_labels))
        self.assertEqual(pdb_has_alt_conf, ampal_has_alt_conf)

    def test_parse_1ek9(self):
        """Check that 3qy1 has been parsed correctly."""
        test_path = os.path.join('unit_tests', 'testing_files', '1ek9.pdb')
        self.check_ampal_contents(test_path)

    def test_parse_2ht0(self):
        """Check that 3qy1 has been parsed correctly."""
        test_path = os.path.join('unit_tests', 'testing_files', '2ht0.pdb')
        self.check_ampal_contents(test_path)

    def test_parse_3qy1(self):
        """Check that 3qy1 has been parsed correctly."""
        test_path = os.path.join('unit_tests', 'testing_files', '3qy1.pdb')
        self.check_ampal_contents(test_path)
