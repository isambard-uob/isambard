import unittest

from hypothesis import given, settings
from hypothesis.strategies import integers, floats, tuples, text
from numpy import allclose

import isambard_dev as isambard


def count_bonds(sequence, phos=False):
    bond_counts = {
        'G': 24,
        'C': 20,
        'A': 23,
        'T': 21
    }
    if not phos:
        return sum([bond_counts[x] for x in sequence]) + (len(sequence) - 1) - 3
    else:
        return sum([bond_counts[x] for x in sequence]) + (len(sequence) - 1)


class TestSingleStrandHelix(unittest.TestCase):

    @given(text('ATGC'))
    @settings(max_examples=50)
    def test_hna_random_seq(self, sequence):
        """Test straight building using random DNA sequences."""
        if not sequence:
            with self.assertRaises(ValueError):
                isambard.ampal.secondary_structure.NucleicAcidStrand(sequence)
        else:
            ssnh = isambard.ampal.secondary_structure.NucleicAcidStrand(sequence)
            self.assertEqual(len(ssnh), len(sequence))

            ideal_bonds = count_bonds(sequence)
            found_bonds = len(isambard.ampal.interactions.find_covalent_bonds(ssnh))
            self.assertEqual(ideal_bonds, found_bonds)

    @given(text('ATGC'))
    @settings(max_examples=50)
    def test_hna_bond_numbers_w_phos(self, sequence):
        if not sequence:
            with self.assertRaises(ValueError):
                isambard.ampal.secondary_structure.NucleicAcidStrand(sequence)
        else:
            ideal_bonds = count_bonds(sequence, phos=True)
            ssnh = isambard.ampal.secondary_structure.NucleicAcidStrand(sequence, phos_3_prime=True)
            found_bonds = len(isambard.ampal.interactions.find_covalent_bonds(ssnh))
            self.assertEqual(ideal_bonds, found_bonds)

    @given(tuples(*[integers(min_value=-10000, max_value=10000) for _ in range(3)]),
           tuples(*[integers(min_value=-10000, max_value=10000) for _ in range(3)]))
    @settings(max_examples=50)
    def test_hna_start_end_int(self, start, end):
        """Test the from_start_and_end class method of the SingleStrandedHelix using ints."""
        sequence = 'GAGATATACACA'
        if start == end:
            with self.assertRaises(ValueError):
                isambard.ampal.secondary_structure.NucleicAcidStrand.from_start_and_end(start, end, sequence)
        else:
            ssnh = isambard.ampal.secondary_structure.NucleicAcidStrand.from_start_and_end(start, end, sequence)
            self.assertEqual(len(ssnh), 12)

    @given(tuples(*[floats(min_value=-10000, max_value=10000) for _ in range(3)]),
           tuples(*[floats(min_value=-10000, max_value=10000) for _ in range(3)]))
    @settings(max_examples=50)
    def test_hna_start_end_floats(self, start, end):
        """Test the from_start_and_end class method of the SingleStrandedHelix using floats."""
        sequence = 'GAGATATACACA'
        if allclose(start, end):
            with self.assertRaises(ValueError):
                isambard.ampal.secondary_structure.NucleicAcidStrand.from_start_and_end(start, end, sequence)
        else:
            ssnh = isambard.ampal.secondary_structure.NucleicAcidStrand.from_start_and_end(start, end, sequence)
            self.assertEqual(len(ssnh), 12)

    @given(text('ATGC'),
           tuples(*[integers(min_value=-10000, max_value=10000) for _ in range(3)]),
           tuples(*[integers(min_value=-10000, max_value=10000) for _ in range(3)]))
    @settings(max_examples=50)
    def test_hna_start_end_ran_seq(self, sequence, start, end):
        """Test SingleStrandHelix with random sequence, start and end."""
        if allclose(start, end) or not sequence:
            with self.assertRaises(ValueError):
                isambard.ampal.secondary_structure.NucleicAcidStrand.from_start_and_end(start, end, sequence)
        else:
            ssnh = isambard.ampal.secondary_structure.NucleicAcidStrand.from_start_and_end(start, end, sequence)
            self.assertEqual(len(ssnh), len(sequence))

            ideal_bonds = count_bonds(sequence)
            found_bonds = len(isambard.ampal.interactions.find_covalent_bonds(ssnh))
            self.assertEqual(ideal_bonds, found_bonds)


class TestDNADuplex(unittest.TestCase):

    @given(text('ATGC'))
    @settings(max_examples=50)
    def test_hna_random_seq(self, sequence):
        """Test straight duplex building using random DNA sequences."""
        if not sequence:
            with self.assertRaises(ValueError):
                isambard.ampal.specifications.DNADuplex.from_sequence(sequence)
        else:
            dd = isambard.ampal.specifications.DNADuplex.from_sequence(sequence)
            self.assertEqual(len(dd), 2)
            self.assertEqual(len(dd[0]), len(sequence))
            self.assertEqual(len(dd[1]), len(sequence))

            ideal_bonds_s1 = count_bonds(sequence)
            ideal_bonds_s2 = count_bonds(
                isambard.ampal.specifications.nucleic_acid_duplex.generate_antisense_sequence(sequence))
            found_bonds_s1 = len(isambard.ampal.interactions.find_covalent_bonds(dd[0]))
            found_bonds_s2 = len(isambard.ampal.interactions.find_covalent_bonds(dd[1]))
            found_bonds_dd = len(isambard.ampal.interactions.find_covalent_bonds(dd))
            self.assertEqual(ideal_bonds_s1, found_bonds_s1)
            self.assertEqual(ideal_bonds_s2, found_bonds_s2)
            self.assertEqual(ideal_bonds_s1 + ideal_bonds_s2, found_bonds_dd)

    @given(tuples(*[integers(min_value=-10000, max_value=10000) for _ in range(3)]),
           tuples(*[integers(min_value=-10000, max_value=10000) for _ in range(3)]))
    @settings(max_examples=50)
    def test_dna_duplex_ints(self, start, end):
        """Test DNADuplex with random int start and end."""
        sequence = 'GAGATATACACA'
        if allclose(start, end):
            with self.assertRaises(ValueError):
                isambard.ampal.specifications.DNADuplex.from_start_and_end(start, end, sequence)
        else:
            dd = isambard.ampal.specifications.DNADuplex.from_start_and_end(start, end, sequence)
            self.assertEqual(len(dd), 2)
            self.assertEqual(len(dd[0]), 12)
            self.assertEqual(len(dd[1]), 12)

    @given(tuples(*[floats(min_value=-10000, max_value=10000) for _ in range(3)]),
           tuples(*[floats(min_value=-10000, max_value=10000) for _ in range(3)]))
    @settings(max_examples=50)
    def test_dna_duplex_floats(self, start, end):
        """Test DNADuplex with random float start and end."""
        sequence = 'GAGATATACACA'
        if allclose(start, end):
            with self.assertRaises(ValueError):
                isambard.ampal.specifications.DNADuplex.from_start_and_end(start, end, sequence)
        else:
            dd = isambard.ampal.specifications.DNADuplex.from_start_and_end(start, end, sequence)
            self.assertEqual(len(dd), 2)
            self.assertEqual(len(dd[0]), 12)
            self.assertEqual(len(dd[1]), 12)

    @given(text('ATGC'),
           tuples(*[integers(min_value=-10000, max_value=10000) for _ in range(3)]),
           tuples(*[integers(min_value=-10000, max_value=10000) for _ in range(3)]))
    @settings(max_examples=50)
    def test_dna_duplex_start_end_ran_seq(self, sequence, start, end):
        """Test SingleStrandHelix with random sequence, start and end."""
        if allclose(start, end) or not sequence:
            with self.assertRaises(ValueError):
                isambard.ampal.specifications.DNADuplex.from_start_and_end(start, end, sequence)
        else:
            dd = isambard.ampal.specifications.DNADuplex.from_start_and_end(start, end, sequence)
            self.assertEqual(len(dd), 2)
            self.assertEqual(len(dd[0]), len(sequence))
            self.assertEqual(len(dd[1]), len(sequence))

            ideal_bonds_s1 = count_bonds(sequence)
            ideal_bonds_s2 = count_bonds(
                isambard.ampal.specifications.nucleic_acid_duplex.generate_antisense_sequence(sequence))
            found_bonds_s1 = len(isambard.ampal.interactions.find_covalent_bonds(dd[0]))
            found_bonds_s2 = len(isambard.ampal.interactions.find_covalent_bonds(dd[1]))
            found_bonds_dd = len(isambard.ampal.interactions.find_covalent_bonds(dd))
            self.assertEqual(ideal_bonds_s1, found_bonds_s1)
            self.assertEqual(ideal_bonds_s2, found_bonds_s2)
            self.assertEqual(ideal_bonds_s1 + ideal_bonds_s2, found_bonds_dd)


__author__ = 'Christopher W. Wood'
