import unittest

import isambard_dev.add_ons.filesystem as filesystem


class NumberOfMmolsTestCase(unittest.TestCase):
    """Tests for isambard.add_ons.filesystem.number_of_mmols"""
    def test_number_of_mmols_4nfw(self):
        """Tests that there are 18 mmol files for the PDB code 4nfw."""
        code = '4nfw'
        self.assertEqual(filesystem.number_of_mmols(code=code), 18)


class PreferredMmolTestCase(unittest.TestCase):
    """Tests isambard.add_ons.filesystem.preferred_mmol """
    def test_preferred_mmol_1cwa(self):
        """ Preferred mmol number for the PDB code 1cwa should be 1 """
        code = '1cwa'
        self.assertEqual(filesystem.preferred_mmol(code=code), 1)

    def test_preferred_mmol_4nfw(self):
        """ Preferred mmol number for the PDB code 1cwa should be 4 """
        code = '4nfw'
        self.assertEqual(filesystem.preferred_mmol(code=code), 4)

    def test_nonsense_raises_value_error(self):
        code = "nonsense_code"
        self.assertRaises(OSError, filesystem.preferred_mmol, code)


class GetMmolTestCase(unittest.TestCase):
    """Tests for isambard.add_ons.filesystem.get_mmol"""
    def test_get_mmol_file_length_4nfw_1(self):
        """The mmol paths in the list should be consistent with calls to get_mmol_path"""
        code = '4nfw'
        self.assertEqual(len(filesystem.get_mmol(code=code, mmol_number=1)), 110328)

    def test_nonsense_is_none(self):
        """A non-existent PDB code should return None. """
        code = "nonsense_code"
        self.assertIsNone(filesystem.get_mmol(code=code))


class PdbeStatusTestCase(unittest.TestCase):
    """Tests for isambard.add_ons.filesystem.pdbe_status_code"""

    def test_status_code_valid(self):
        # Test on a valid pdb code, e.g. 1cwa.
        code = '1cwa'
        self.assertEqual(filesystem.pdbe_status_code(code), 200)

    def test_status_code_invalid(self):
        code = 'invalid_code'
        self.assertEqual(filesystem.pdbe_status_code(code), 404)


__author__ = 'Jack W. Heal'
