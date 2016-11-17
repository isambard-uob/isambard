import unittest

from hypothesis import given, settings
from hypothesis.strategies import text

import isambard_dev as isambard


class TestScwrl4(unittest.TestCase):

    @given(text('ACDEFGHIKLMNPQRSTVWY', min_size=5, max_size=5))
    @settings(max_examples=20)
    def test_helix_seq_packing(self, sequence):
        """Test Scwrl sidechain packing using a small helix."""
        helix = isambard.specifications.Helix(5)
        helix.pack_new_sequence(sequence)
        self.assertEqual(helix.sequence, sequence)

    @given(text('ABCDEFGHIJKLMNOPQRSTUVWXYZ', min_size=5, max_size=5))
    @settings(max_examples=20)
    def test_helix_non_canonical(self, sequence):
        """Test Scwrl sidechain packing using a small helix."""
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
