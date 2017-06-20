import glob
import os
import unittest

import isambard_dev as isambard
import isambard_dev.add_ons.pacc as pacc


class PACCTestCase(unittest.TestCase):
    """Tests for isambard.add_ons.pacc"""

    def setUp(self):
        test_cc_folder = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', 'test_ccs')
        self.pdb_paths = glob.glob('{}/*'.format(test_cc_folder))

    def check_register(self, heptad_register, sequence_len):
        self.assertEqual(len(heptad_register), sequence_len)
        return

    def check_parameters(self, parameter_values, sequence_len):
        self.assertEqual(len(parameter_values), sequence_len)
        return

    def test_pacc_on_ccs(self):
        """Tests PACC for each file in testing_files/test_ccs"""
        for file in self.pdb_paths:
            cc = isambard.ampal.convert_pdb_to_ampal(file)
            ccr = isambard.ampal.Assembly()
            for h in cc:
                ccr.append(h[1:28])
            cc_pacc = pacc.PACCAnalysis(ccr)
            # Check register assignment
            register, fit = cc_pacc.heptad_register()
            self.check_register(register, cc_pacc.cc_len)
            for p in [cc_pacc.radii_layers, cc_pacc.alpha_layers, cc_pacc.ca_layers]:
                self.check_parameters(cc_pacc.calc_average_parameters(p)[0], cc_pacc.cc_len)
        return
