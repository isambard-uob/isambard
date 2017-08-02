import copy
import unittest

import isambard_dev as isambard
import numpy


@unittest.skipUnless(
    'scwrl' in isambard.settings.global_settings,
    "External program not detected.")
class TestConvertProToHyp(unittest.TestCase):
    def test_convert_pro_to_hyp(self):
        col = isambard.specifications.CoiledCoil.tropocollagen(aa=21)
        col.pack_new_sequences(['GPPGPPGPPGPPGPPGPPGPP'] * 3)
        to_convert = [
            res for (i, res) in enumerate(col.get_monomers())
            if not (i + 1) % 3]
        ori_pros = copy.deepcopy(to_convert)
        for pro in to_convert:
            isambard.ampal.non_canonical.convert_pro_to_hyp(pro)
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
