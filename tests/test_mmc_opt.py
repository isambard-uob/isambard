"""Tests for MMC Optimiser."""

import unittest

from hypothesis import given
from hypothesis.strategies import floats
import isambard
import isambard.specifications
from isambard.optimisation.mmc_optimizer import (
    MMCParameter, MMCParameterType, MMCParameterOptimisation)


class TestMMCParameter(unittest.TestCase):
    """Tests parameter creation and randomisation."""
    @given(floats(allow_nan=False, allow_infinity=False),
           floats(allow_nan=False, allow_infinity=False))
    def test_unifrom_parameter(self, val_a, val_b):
        """Make new uniform parameter and checks value."""
        uniform = MMCParameter('Uniform', MMCParameterType.UNIFORM_DIST,
                               (val_a, val_b))
        uniform.randomise_proposed_value()
        if ((val_a - val_b) == float('inf')) or (
                (val_a - val_b) == float('-inf')):
            pass
        elif val_a < val_b:
            self.assertTrue(val_a <= uniform.proposed_value <= val_b)
        else:
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

    def test_build_from_parameters(self):
        """Tests model building using parameters."""
        input_parameters = (
            MMCParameter('Oligomeric State', MMCParameterType.STATIC_VALUE, 2),
            MMCParameter('Amino Acids', MMCParameterType.STATIC_VALUE, 28),
            MMCParameter('Radius', MMCParameterType.UNIFORM_DIST, (2, 10)),
            MMCParameter('Pitch', MMCParameterType.UNIFORM_DIST, (50, 300)),
            MMCParameter('PhiCA', MMCParameterType.UNIFORM_DIST, (150, 360)),
        )
        build_args = [p.current_value for p in input_parameters]
        coco = isambard.specifications.CoiledCoil.from_parameters(*build_args)
        self.assertEqual(len(coco), 2)
        self.assertEqual([len(x) for x in coco], [28, 28])
        self.assertTrue(all([2 <= rad <= 10 for rad in coco.major_radii]))
        self.assertTrue(all([50 <= rad <= 300 for rad in coco.major_pitches]))
        self.assertTrue(all([150 <= rad <= 360 for rad in coco.phi_c_alphas]))


class TestMMCOptimizer(unittest.TestCase):
    """Tests functionality in the MMC optimizer."""
    @given(floats(allow_nan=False, allow_infinity=False),
           floats(allow_nan=False, allow_infinity=False),
           floats(allow_nan=False, allow_infinity=False))
    def test_check_move(self, old, new, temp):
        """Checks a random move."""
        MMCParameterOptimisation.check_move(old, new, temp)

    @given(floats(allow_nan=False, allow_infinity=False),
           floats(allow_nan=False, allow_infinity=False))
    def test_check_move_frozen(self, old, new):
        """Checks a move at 0 K, which should always be false."""
        self.assertFalse(MMCParameterOptimisation.check_move(old, new, 0))
