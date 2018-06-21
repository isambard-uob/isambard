import glob
import os
import unittest

import isambard
import ampal
import isambard.evaluation.pacc as pacc
import numpy as np

class PACCTestCase(unittest.TestCase):
    """Tests for isambard.evaluation.pacc"""

    def setUp(self):
        self.test_cc_folder = os.path.join(os.path.dirname(__file__), 'testing_files', 'test_pacc')
        self.pdb_paths = glob.glob('{}/*'.format(self.test_cc_folder))

    def check_register(self, heptad_register, sequence_len):
        self.assertEqual(len(heptad_register), sequence_len)
        return

    def check_parameters(self, parameter_values, sequence_len):
        self.assertEqual(len(parameter_values), sequence_len)
        return

    def test_pacc_on_ccs(self):
        """Tests PACC for each file in testing_files/test_ccs"""
        assert len(self.pdb_paths) > 0
        for file in self.pdb_paths:
            cc = ampal.load_pdb(file)
            ccr = ampal.Assembly()
            for h in cc:
                ccr.append(h[1:28])
            cc_pacc = pacc.PACCAnalysis(ccr)
            # Check register assignment
            register, fit = cc_pacc.heptad_register()
            self.check_register(register, cc_pacc.cc_len)
            for p in [cc_pacc.radii_layers, cc_pacc.alpha_layers, cc_pacc.ca_layers]:
                self.check_parameters(cc_pacc.calc_average_parameters(p)[0], cc_pacc.cc_len)
        return

    def test_pacc_parallel_antiparallel(self):
        file_name = os.path.join(self.test_cc_folder, 'p3_p4.pdb')
        cc = ampal.load_pdb(file_name)
        cca = pacc.PACCAnalysis(cc)
        radius = np.mean(cca.radii_layers)
        assert radius < 5.1

        file_name = os.path.join(self.test_cc_folder, 'APH.pdb')
        cc = ampal.load_pdb(file_name)
        cca = pacc.PACCAnalysis(cc)
        radius = np.mean(cca.radii_layers)
        assert radius < 5.1
