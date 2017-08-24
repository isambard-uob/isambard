"""Tests for MMC Optimiser."""

import unittest

from hypothesis import given
from hypothesis.strategies import floats
from isambard_dev.optimisation.mmc_optimizer import (
    MMCParameter, MMCParameterType)


class TestMMCParameter(unittest.TestCase):
    """Tests parameter creation and randomisation."""
    @given(floats(allow_nan=False, allow_infinity=False),
           floats(allow_nan=False, allow_infinity=False))
    def test_unifrom_parameter(self, val_a, val_b):
        """Make new uniform parameter and checks value."""
        uniform = MMCParameter('Uniform', MMCParameterType.UNIFORM_DIST,
                               (val_a, val_b))
        uniform.randomise_proposed_value()
        if val_a < val_b:
            self.assertTrue(val_a <= uniform.proposed_value <= val_b)
        if val_a > val_b:
            self.assertTrue(val_a >= uniform.proposed_value >= val_b)

    @given(floats(allow_nan=False, allow_infinity=False),
           floats(allow_nan=False, allow_infinity=False))
    def test_normal_parameter(self, val_a, val_b):
        """Make new normal parameter and checks value."""
        normal = MMCParameter('Normal', MMCParameterType.NORMAL_DIST,
                              (val_a, val_b))
        normal.randomise_proposed_value()
        self.assertTrue(normal.proposed_value is not None)

    def test_accept_parameter(self):
        """Tests parameter acceptance."""
        uniform = MMCParameter('Uniform', MMCParameterType.UNIFORM_DIST,
                               (100, 60))
        self.assertTrue(uniform.current_value is not None)
        self.assertTrue(uniform.proposed_value is None)
        uniform.randomise_proposed_value()
        prop_val = uniform.proposed_value
        self.assertTrue(uniform.proposed_value is not None)
        uniform.accept_proposed_value()
        self.assertTrue(uniform.current_value == prop_val)
        self.assertTrue(uniform.proposed_value is None)

    def test_reject_parameter(self):
        """Tests parameter rejection."""
        uniform = MMCParameter('Uniform', MMCParameterType.UNIFORM_DIST,
                               (100, 60))
        cur_val = uniform.current_value
        self.assertTrue(uniform.current_value is not None)
        self.assertTrue(uniform.proposed_value is None)
        uniform.randomise_proposed_value()
        prop_val = uniform.proposed_value
        self.assertTrue(uniform.current_value == cur_val)
        self.assertTrue(uniform.proposed_value is not None)
        uniform.reject_proposed_value()
        self.assertTrue(uniform.current_value == cur_val)
        self.assertTrue(uniform.current_value != prop_val)
        self.assertTrue(uniform.proposed_value is None)
