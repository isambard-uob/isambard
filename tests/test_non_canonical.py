"""Tests functionality in ampal/non_cannonical.py"""

import copy
import unittest

import isambard
import isambard.specifications as specs
from isambard.modelling import convert_pro_to_hyp
from isambard.modelling.scwrl import scwrl_available, pack_side_chains_scwrl
import numpy


@unittest.skipUnless(scwrl_available(), "External program not detected.")
class TestConvertProToHyp(unittest.TestCase):
    def test_convert_pro_to_hyp(self):
        col = specs.CoiledCoil.tropocollagen(aa=21)
        col = pack_side_chains_scwrl(col, ['GPPGPPGPPGPPGPPGPPGPP'] * 3)
        to_convert = [
            res for (i, res) in enumerate(col.get_monomers())
            if not (i + 1) % 3]
        ori_pros = copy.deepcopy(to_convert)
        for pro in to_convert:
            convert_pro_to_hyp(pro)
        self.assertEqual(col.sequences, ['GPXGPXGPXGPXGPXGPXGPX'] * 3)
        hyps = to_convert
        self.assertEqual(len(hyps), 7 * 3)
        hyp_atom_labels = ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OD1')
        common_atoms = ('N', 'CA', 'C', 'O', 'CB')
        for (pro, hyp) in zip(ori_pros, hyps):
            for (label, atom) in hyp.atoms.items():
                self.assertTrue(label in hyp_atom_labels)
                if label in common_atoms:
                    self.assertTrue(
                        numpy.allclose(atom.array, hyp[label].array))
