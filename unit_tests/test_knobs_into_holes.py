import unittest
import os

import isambard
import isambard.add_ons.knobs_into_holes as knobs_into_holes
import isambard.add_ons.parmed_to_ampal as parmed_to_ampal


class KnobGroupTestCase(unittest.TestCase):
    """Tests for class KnobGroup"""
    def setUp(self):
        test_file = os.path.join('unit_tests', 'testing_files', '2ebo_1.mmol')
        self.test_assembly = isambard.ampal.convert_pdb_to_ampal(test_file)
        test_file = os.path.join('unit_tests', 'testing_files', '2j58_1.cif')
        self.test_assembly_2 = parmed_to_ampal.convert_cif_to_ampal(test_file)

    def test_number_of_kihs(self):
        """ Test there are 38 kihs found at 7.0 cutoff for 2ebo. """
        a = self.test_assembly
        kg = knobs_into_holes.KnobGroup.from_helices(a, min_helix_length=0)
        self.assertTrue(len(kg) == 38)
        kihs_locations = ''
        for i in range(len(a.helices)):
            kihs_locations += str(len([x for x in kg if x.knob_helix.number == i]))
        self.assertTrue(int(kihs_locations) == 914804714)

    def test_kihs_at_large_cutoffs(self):
        for scut in [7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]:
            a = self.test_assembly
            kg = knobs_into_holes.KnobGroup.from_helices(a, cutoff=scut)
            self.assertEqual(len(kg), len(set([str(x) for x in kg])))