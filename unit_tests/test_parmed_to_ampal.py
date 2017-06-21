import os
import unittest

import parmed

import isambard_dev as isambard
import isambard_dev.add_ons.parmed_to_ampal as parmed_to_ampal


class ParmedToAmpalTestCase(unittest.TestCase):
    """ tests for isambard.add_ons.parmed_to_ampal"""
    def setUp(self):
        test_files = [os.path.join('unit_tests', 'testing_files', x)
                      for x in ['1ek9.pdb', '2ht0.pdb', '3qy1.pdb']]
        self.test_ampals = [isambard.ampal.convert_pdb_to_ampal(x) for x in test_files]
        test_parmeds = [parmed.formats.PDBFile.parse(x, skip_bonds=True) for x in test_files]
        self.test_parmeds = [parmed_to_ampal.parmed_to_ampal(x) for x in test_parmeds]

    def test_number_of_atoms(self):
        tests = []
        for ta, tp in zip(self.test_ampals, self.test_parmeds):
            at = len(list(ta.get_atoms())) == len(list(tp.get_atoms()))
            tests.append(at)
        self.assertTrue(all(tests))

    def test_number_of_residues(self):
        tests = []
        for ta, tp in zip(self.test_ampals, self.test_parmeds):
            at = len(list(ta.get_monomers())) == len(list(tp.get_monomers()))
            tests.append(at)
        self.assertTrue(all(tests))

    def test_length_of_pdb_attribute(self):
        tests = []
        for ta, tp in zip(self.test_ampals, self.test_parmeds):
            at = len(ta.pdb) == len(tp.pdb)
            tests.append(at)
        self.assertTrue(all(tests))

    def test_polymer_types(self):
        tests = []
        for ta, tp in zip(self.test_ampals, self.test_parmeds):
            at1 = len(ta) == len(tp)
            at2 = len([x for x in ta if isinstance(x, isambard.ampal.Polypeptide)]) == \
                  len([x for x in tp if isinstance(x, isambard.ampal.Polypeptide)])
            at3 = len([x for x in ta if isinstance(x, isambard.ampal.Polynucleotide)]) == \
                  len([x for x in tp if isinstance(x, isambard.ampal.Polynucleotide)])
            tests.append(all([at1, at2, at3]))
        self.assertTrue(all(tests))


__author__ = 'Jack W. Heal'
