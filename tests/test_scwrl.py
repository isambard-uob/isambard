import os
import unittest

import ampal
from hypothesis import given, settings
from hypothesis.strategies import integers, text

import isambard.specifications as specs
from isambard.modelling.scwrl import scwrl_available, pack_side_chains_scwrl


cannonical_labels = 'ACDEFGHIKLMNPQRSTVWY'
non_cannonical_labels = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'


@unittest.skipUnless(scwrl_available(), "External program not detected.")
class TestScwrl4(unittest.TestCase):

    @given(text(cannonical_labels, min_size=5, max_size=5))
    @settings(max_examples=20)
    def test_helix_seq_packing(self, sequence):
        """Test Scwrl sidechain packing using a small helix."""
        helix = ampal.Assembly(specs.Helix(5))
        helix = pack_side_chains_scwrl(helix, [sequence])
        self.assertEqual(helix.sequences, [sequence])

    @given(text(non_cannonical_labels, min_size=5, max_size=5))
    @settings(max_examples=20)
    def test_helix_non_canonical(self, sequence):
        """Test Scwrl sidechain packing when given non-cannonical residue labels."""
        helix = ampal.Assembly(specs.Helix(5))
        helix = pack_side_chains_scwrl(helix, [sequence])
        non_can = [
            x for x in sequence if x not in cannonical_labels]
        r_seq = sequence
        if non_can:
            for aa in non_can:
                r_seq = r_seq.replace(aa, 'G')
            self.assertEqual(helix.sequences, [r_seq])
        else:
            self.assertEqual(helix.sequences, [r_seq])

    @given(text(non_cannonical_labels + non_cannonical_labels.lower(), min_size=5, max_size=5))
    @settings(max_examples=20)
    def test_helix_non_canonical_lower(self, sequence):
        """Test Scwrl sidechain packing using upper and lower case labels."""
        helix = ampal.Assembly(specs.Helix(5))
        helix = pack_side_chains_scwrl(helix, [sequence])
        r_seq = sequence[:]
        non_can = [x for x in sequence if x.upper() not in cannonical_labels]
        if non_can:
            for aa in non_can:
                r_seq = r_seq.replace(aa, 'G')
            self.assertEqual(helix.sequences, [r_seq.upper()])
        else:
            self.assertEqual(helix.sequences, [r_seq.upper()])

    def test_basis_set_dimer(self):
        """Test packing basis set dimer model."""
        cc_di = specs.CoiledCoil(2)
        cc_packed = pack_side_chains_scwrl(cc_di, cc_di.basis_set_sequences)
        self.assertEqual(cc_packed.sequences, cc_di.basis_set_sequences)

    def test_basis_set_trimer(self):
        """Test packing basis set trimer model."""
        cc_tri = specs.CoiledCoil(3)
        cc_packed = pack_side_chains_scwrl(cc_tri, cc_tri.basis_set_sequences)
        self.assertEqual(cc_packed.sequences, cc_tri.basis_set_sequences)

    def test_basis_set_tetramer(self):
        """Test packing basis set tetramer model."""
        cc_tet = specs.CoiledCoil(4)
        cc_packed = pack_side_chains_scwrl(cc_tet, cc_tet.basis_set_sequences)
        self.assertEqual(cc_packed.sequences, cc_tet.basis_set_sequences)
