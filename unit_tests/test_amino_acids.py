import unittest

import isambard


class IdealBackboneBondAnglesTestCase(unittest.TestCase):
    """ Tests backbone_bond_angles dictionary"""

    def test_c_angles(self):
        """ Test the angles around the backbone C atom add to 360. """
        trans = isambard.tools.amino_acids.ideal_backbone_bond_angles['trans']
        cis = isambard.tools.amino_acids.ideal_backbone_bond_angles['cis']
        sum_1 = trans['ca_c_n'] + trans['ca_c_o'] + trans['o_c_n']
        sum_2 = cis['ca_c_n'] + cis['ca_c_o'] + cis['o_c_n']
        self.assertEqual(sum_1, sum_2)
        self.assertEqual(sum_1, 360.0)


class StandardAminoAcidsTestCase(unittest.TestCase):
    """Tests isambard.tools.amino_acids.tools_standard_amino_acids """

    def test_key_length(self):
        """All standard amino acid codes should be strings of length 3"""
        key_lengths = [len(k) for k in isambard.tools.amino_acids.standard_amino_acids.keys()]
        self.assertTrue(all(x == 1 for x in key_lengths))

    def test_value_length(self):
        """All standard amino acid letters should be strings of length 1"""
        val_lengths = [len(v) for v in isambard.tools.amino_acids.standard_amino_acids.values()]
        self.assertTrue(all(x == 3 for x in val_lengths))


class AACodeTestCase(unittest.TestCase):
    """Tests for isambard.tools.amino_acids.tools_get_aa_code"""

    def test_aa_code_invalid(self):
        self.assertIsNone(isambard.tools.amino_acids.get_aa_code('testingnonsensestring'))
        self.assertIsNone(isambard.tools.amino_acids.get_aa_code('X'))

    def test_aa_code_valid(self):
        self.assertEqual(isambard.tools.amino_acids.get_aa_code('A'), 'ALA')
        self.assertEqual(isambard.tools.amino_acids.get_aa_code('C'), 'CYS')


class AALetterTestCase(unittest.TestCase):
    """Tests for isambard.tools.amino_acids.tools_get_aa_letter"""

    def test_aa_letter_invalid(self):
        self.assertEqual(isambard.tools.amino_acids.get_aa_letter('longstringofletters'), 'X')

    def test_aa_letter_valid(self):
        self.assertEqual(isambard.tools.amino_acids.get_aa_letter('ALA'), 'A')
        self.assertEqual(isambard.tools.amino_acids.get_aa_letter('TYR'), 'Y')


class AminoAcidDictTestCase(unittest.TestCase):
    """Tests isambard.tools.amino_acids.tools_amino_acid_dict"""

    def test_standard_amino_acids(self):
        """Check that there are 20 standard amino acids in amino_acids_dict."""
        standards = [x['letter'] for x in isambard.tools.amino_acids.amino_acids_dict.values() if x['letter'] != "X"]
        self.assertEqual(len(standards), 20)

    def test_values_type(self):
        """Check all values in amino_acids_dict are dictionaries."""
        self.assertTrue(all(isinstance(x, dict) for x in isambard.tools.amino_acids.amino_acids_dict.values()))

    def test_unexpected_letters(self):
        unexpected_letters = ['B', 'J', 'O', 'U', 'Z']
        self.assertEqual(len(
            [x for x in isambard.tools.amino_acids.amino_acids_dict.values() if x['letter'] in unexpected_letters]), 0)


class GetAAInfoTestCase(unittest.TestCase):
    """Tests isambard.tools.amino_acids.tools_get_aa_info"""

    def test_nonsense_code(self):
        self.assertRaises(IOError, isambard.tools.amino_acids.get_aa_info, code="nonsensecode")

    def test_phi(self):
        aad = isambard.tools.amino_acids.get_aa_info("PHI")
        self.assertIs(type(aad), dict)
        self.assertEqual(aad['description'], 'IODO-PHENYLALANINE')
        self.assertEqual(aad['modified'], 'PHE')


class AddAminoAcidToJsonTestCase(unittest.TestCase):
    """Tests isambard.tools.amino_acids.tools_add_amino_acid_to_json"""

    def test_exisitng_amino_acid(self):
        self.assertRaises(IOError, isambard.tools.amino_acids.add_amino_acid_to_json,
                          code='PHE', description='Phenylalanine')


__author__ = 'Jack W. Heal'
