"""Tests for the dssp.py module."""

import unittest
import warnings

from isambard.evaluation.dssp import dssp_available, tag_dssp_data
from isambard.modelling.scwrl import pack_side_chains_scwrl
import isambard.specifications as specs

warnings.filterwarnings("ignore")


@unittest.skipUnless(dssp_available(), "External program not detected.")
class TestDSSP(unittest.TestCase):
    """Tests the interface to DSSP."""

    def test_dimer_tags(self):
        """Checks the DSSP data tagged to a coiled-coil dimer."""
        cc_di = specs.CoiledCoil(2)
        cc_packed = pack_side_chains_scwrl(cc_di, cc_di.basis_set_sequences)
        tag_dssp_data(cc_packed)
        for residue in cc_packed.get_monomers():
            self.assertTrue('dssp_data' in residue.tags)
