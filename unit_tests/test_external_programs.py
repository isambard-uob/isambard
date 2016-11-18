import os
import unittest
import warnings

from hypothesis import given, settings
from hypothesis.strategies import integers, text

import isambard_dev as isambard

warnings.filterwarnings("ignore")


cannonical_labels = 'ACDEFGHIKLMNPQRSTVWY'
non_cannonical_labels = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'


class TestScwrl4(unittest.TestCase):

    def test_check_scwrl_avail(self):
        """Test if ISAMBARD can detect Scwrl, check its availability and behave accordingly."""
        old_path = isambard.settings.global_settings['scwrl']['path']
        try:
            avail = isambard.external_programs.scwrl.check_scwrl_avail()
            self.assertTrue(avail)

            helix = isambard.specifications.Helix(30)
            scwrl_out = isambard.external_programs.scwrl.pack_sidechains(helix.pdb, 'V'*30, path=False)
            self.assertEqual(type(scwrl_out[0]), str)
            self.assertEqual(type(scwrl_out[1]), float)

            isambard.settings.global_settings['scwrl']['path'] = ''
            avail = isambard.external_programs.scwrl.check_scwrl_avail()
            isambard.settings.global_settings['scwrl']['available'] = avail
            self.assertFalse(avail)

            scwrl_out = isambard.external_programs.scwrl.pack_sidechains(helix.pdb, 'V'*30, path=False)
            self.assertIsNone(scwrl_out)
        finally:
            isambard.settings.global_settings['scwrl']['path'] = old_path
            isambard.settings.global_settings[
                'scwrl']['available'] = isambard.external_programs.scwrl.check_scwrl_avail()

    @given(text(cannonical_labels, min_size=5, max_size=5))
    @settings(max_examples=20)
    def test_helix_seq_packing(self, sequence):
        """Test Scwrl sidechain packing using a small helix."""
        helix = isambard.specifications.Helix(5)
        helix.pack_new_sequence(sequence)
        self.assertEqual(helix.sequence, sequence)

    @given(text(non_cannonical_labels, min_size=5, max_size=5))
    @settings(max_examples=20)
    def test_helix_non_canonical(self, sequence):
        """Test Scwrl sidechain packing when given non-cannonical residue labels."""
        helix = isambard.specifications.Helix(5)
        helix.pack_new_sequence(sequence)
        r_seq = sequence[:]
        non_can = [x for x in sequence if x not in isambard.tools.amino_acids.standard_amino_acids.keys()]
        if non_can:
            for aa in non_can:
                r_seq = r_seq.replace(aa, 'G')
            self.assertEqual(helix.sequence, r_seq)
        else:
            self.assertEqual(helix.sequence, r_seq)

    @given(text(non_cannonical_labels + non_cannonical_labels.lower(), min_size=5, max_size=5))
    @settings(max_examples=20)
    def test_helix_non_canonical_lower(self, sequence):
        """Test Scwrl sidechain packing using upper and lower case labels."""
        helix = isambard.specifications.Helix(5)
        helix.pack_new_sequence(sequence)
        r_seq = sequence[:]
        non_can = [x for x in sequence if x.upper() not in isambard.tools.amino_acids.standard_amino_acids.keys()]
        if non_can:
            for aa in non_can:
                r_seq = r_seq.replace(aa, 'G')
            self.assertEqual(helix.sequence, r_seq.upper())
        else:
            self.assertEqual(helix.sequence, r_seq.upper())

    def test_basis_set_dimer(self):
        """Test packing basis set dimer model."""
        cc_di = isambard.specifications.CoiledCoil(2)
        cc_di.pack_new_sequences(cc_di.basis_set_sequences)
        self.assertEqual(cc_di.sequences, cc_di.basis_set_sequences)

    def test_basis_set_trimer(self):
        """Test packing basis set trimer model."""
        cc_tri = isambard.specifications.CoiledCoil(3)
        cc_tri.pack_new_sequences(cc_tri.basis_set_sequences)
        self.assertEqual(cc_tri.sequences, cc_tri.basis_set_sequences)

    def test_basis_set_tetramer(self):
        """Test packing basis set tetramer model."""
        cc_tet = isambard.specifications.CoiledCoil(4)
        cc_tet.pack_new_sequences(cc_tet.basis_set_sequences)
        self.assertEqual(cc_tet.sequences, cc_tet.basis_set_sequences)


class TestDSSP(unittest.TestCase):

    def check_dssp_tag(self, test_file_path):
        ampal = isambard.ampal.convert_pdb_to_ampal(test_file_path)
        ampal.tag_secondary_structure()
        ss_log = []
        for mon in ampal.get_monomers(ligands=False):
            ss_log.append('secondary_structure' in mon.tags)
        print(ss_log)
        self.assertTrue(all(ss_log))

    def test_check_dssp_avail(self):
        """Test if ISAMBARD can detect DSSP, check its availability and behave accordingly."""
        old_path = isambard.settings.global_settings['dssp']['path']
        try:
            avail = isambard.external_programs.dssp.check_dssp_avail()
            self.assertTrue(avail)
            helix = isambard.specifications.Helix(30)

            dssp_out_ampal = isambard.external_programs.dssp.run_dssp(helix.pdb, path=False)
            self.assertEqual(type(dssp_out_ampal), str)

            isambard.settings.global_settings['dssp']['path'] = ''
            avail = isambard.external_programs.dssp.check_dssp_avail()
            isambard.settings.global_settings['dssp']['available'] = avail
            self.assertFalse(avail)

            dssp_out_ampal = isambard.external_programs.dssp.run_dssp(helix.pdb, path=False)
            self.assertIsNone(dssp_out_ampal)
        finally:
            isambard.settings.global_settings['dssp']['path'] = old_path
            isambard.settings.global_settings['dssp']['available'] = isambard.external_programs.dssp.check_dssp_avail()

    @given(integers(min_value=6, max_value=100))
    @settings(max_examples=20)
    def test_helix_dssp(self, hel_len):
        """Test DSSP helix finding ability."""
        helix = isambard.specifications.Helix(hel_len)
        helical = helix.helices
        self.assertEqual(len(helical), 1)
        self.assertEqual(len(list(helical.get_monomers())), hel_len - 2)

    def test_tag_secondary_structure_3qyi(self):
        """Test the SS tagging functionality."""
        test_file_path = os.path.join('unit_tests', 'testing_files', '3qy1.pdb')
        self.check_dssp_tag(test_file_path)
