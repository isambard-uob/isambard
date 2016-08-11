import random
import unittest
import numpy

from hypothesis import given, settings
from hypothesis.strategies import floats, tuples, lists

import isambard_dev as isambard

from unit_tests.random_isambard_objects import random_vectors, random_floats
from unit_tests.random_isambard_objects import random_angles


class DihedralTestCase(unittest.TestCase):
    """Tests for isambard.tools.dihedral"""

    def test_dihedral_90(self):
        test_points_90 = [(0, 0, 1), (0, 0, 0), (0, 1, 0), (1, 1, 0)]
        self.assertEqual(isambard.geometry.dihedral(*test_points_90), 90.0)
        test_points_90 = [(0, 0, -1), (0, 0, 0), (0, -1, 0), (-1, -1, 0)]
        self.assertEqual(isambard.geometry.dihedral(*test_points_90), -90.0)

    def test_dihedral_0(self):
        test_points_0 = [(0, 0, 1), (0, 0, 0), (0, 1, 0), (0, 1, 1)]
        self.assertEqual(isambard.geometry.dihedral(*test_points_0), 0.0)

    def test_dihedral_180(self):
        test_points_180 = [(0, 0, 1), (0, 0, 0), (0, 1, 0), (0, 1, -1)]
        self.assertEqual(abs(isambard.geometry.dihedral(*test_points_180)), 180.0)

    def test_dihedral_floats(self):
        test_points = (12.1, 23.5, 3.3), (11.1, 24.5, 2.3), (8.1, 3.5, 0.2), (19.2, -5.5, -6.3)
        self.assertAlmostEqual(isambard.geometry.dihedral(*test_points), -68.74, places=2)

    def test_dihedral_same_point(self):
        test_points = [(12.1, 23.5, 3.3)] * 4
        self.assertEqual(isambard.geometry.dihedral(*test_points), 0.0)

    def test_dihedral_range(self):
        """Test dihedral angle is in range [-180.0, 180.0] for randomly generated sets of four points. Repeat x 1000."""
        results = []
        for _ in range(1000):
            # Generate four random points in an array.
            points = numpy.random.rand(4, 3)
            # Multiply them randomly by some integer in the range [-100, 100].
            choices = numpy.random.choice(range(-100, 101), size=12).reshape(4, 3)
            points = numpy.multiply(points, choices)
            # Calculate dihedral angle using the random points and check it lies in correct range
            angle = isambard.geometry.dihedral(*points)
            if (angle <= 180.0) and (angle >= -180.0):
                results.append(1)
            else:
                results.append(0)

        at = numpy.all(results)
        self.assertTrue(at)
        return

    @given(tuples(lists(floats(allow_nan=False, allow_infinity=False), min_size=3, max_size=3),
                  lists(floats(allow_nan=False, allow_infinity=False), min_size=3, max_size=3),
                  lists(floats(allow_nan=False, allow_infinity=False), min_size=3, max_size=3),
                  lists(floats(allow_nan=False, allow_infinity=False), min_size=3, max_size=3)))
    @settings(max_examples=50)
    def test_dihedral_reverves_vectors(self, vectors):
        a, b, c, d = vectors
        angle_1 = isambard.geometry.dihedral(a, b, c, d)
        angle_2 = isambard.geometry.dihedral(d, c, b, a)
        at = numpy.isclose(angle_1, angle_2) or numpy.isclose(angle_1, -angle_2)
        self.assertTrue(at)


class DistanceTestCase(unittest.TestCase):
    """Tests for isambard.tools.geometry.distance"""

    def test_distance(self):
        for z in range(-100, 105, 5):
            self.assertEqual(isambard.geometry.distance((0, 0, 0), (0, 0, z)), abs(z))

    def test_distance_float(self):
        self.assertAlmostEqual(isambard.geometry.distance(
            (11.43, -12.9, -23.11), (-32.43, 42.9, 33.11)), 90.54, places=2)

    def test_distance_in_unit_cube(self):
        """ Generate pair of random points in unit cube, check distance lies within bounds. Repeat x 1000."""
        results = []
        # maximum distance between two points on unit cube is sqrt(3) (e.g. [0, 0, 1] and [1, 1, 0]
        max_d = numpy.sqrt(3)
        for _ in range(1000):
            points = numpy.random.rand(2, 3)
            d = isambard.geometry.distance(*points)
            if (d >= 0.0) and (d <= max_d):
                results.append(1)
            else:
                results.append(0)
        at = numpy.all(results)
        self.assertTrue(at)
        return

    @given(tuples(lists(floats(allow_nan=False, allow_infinity=False), min_size=3, max_size=3),
                  lists(floats(allow_nan=False, allow_infinity=False), min_size=3, max_size=3)))
    @settings(max_examples=50)
    def test_distance_reflexive(self, vectors):
        a, b = vectors
        d1 = isambard.geometry.distance(a, b)
        d2 = isambard.geometry.distance(b, a)
        self.assertEqual(d1, d2)


class UnitVectorTestCase(unittest.TestCase):
    """Tests for isambard.tools.geometry.unit_vector"""

    def test_origin_error(self):
        self.assertRaises(ZeroDivisionError, isambard.geometry.unit_vector, [0, 0, 0])

    @given(
        lists(floats(min_value=1e-100, max_value=1e100, allow_nan=False, allow_infinity=False),
              min_size=3, max_size=3))
    @settings(max_examples=1000)
    def test_unit_length(self, vector):
        """Check norm of the unit_vector for a random array is 1.0. Repeat x 1000."""
        uv = isambard.geometry.unit_vector(vector)
        norm = numpy.linalg.norm(uv)
        self.assertTrue(numpy.isclose(norm, 1.0))


class AngleBetweenVectorsTestCase(unittest.TestCase):
    """Tests for isambard.tools.geometry.angle_between_vectors"""
    @given(tuples(
        lists(floats(min_value=1e-10, max_value=1e-10, allow_nan=False, allow_infinity=False), min_size=3,
              max_size=3),
        lists(floats(min_value=1e-10, max_value=1e-10, allow_nan=False, allow_infinity=False), min_size=3,
              max_size=3)))
    @settings(max_examples=1000)
    def test_angle_range(self, vectors):
        """Test return value is in range [0.0, 180.0] for randomly generated pairs of vectors. Repeat x 1000."""
        a, b = vectors
        angle = isambard.geometry.angle_between_vectors(a, b)
        self.assertTrue((angle >= 0) and (angle <= 180))


class RadiusCCTestCase(unittest.TestCase):
    """Tests for isambard.tools.tool_geometry.radius_of_circumcircle"""

    def test_radius_of_circumcircle(self):
        self.assertAlmostEqual(isambard.geometry.radius_of_circumcircle((-1, 0, 0), (0, 1, 0), (1, 0, 0)), 1.0)
        return

    def test_radius_of_circumcircle_floats(self):
        self.assertAlmostEqual(isambard.geometry.radius_of_circumcircle((-1.5, 0, 0), (0, 1.5, 0), (1.5, 0, 0)), 1.5)
        return


class RotationMatrixTestCase(unittest.TestCase):
    """Tests for isambard.tools.tool_geometry.rotation_matrix"""

    def test_rotation_90(self):
        r_matrix_90 = isambard.geometry.rotation_matrix(numpy.pi / 2, numpy.array([0, 0, 1]), numpy.array([0, 0, 0]))
        aet = numpy.allclose(numpy.dot(r_matrix_90, numpy.array([1, 0, 0, 0])), numpy.array([0, 1, 0, 0]))
        self.assertTrue(aet)


class ClosestDistanceTestCase(unittest.TestCase):
    """Tests for isambard.tools.tool_geometry.closest_distance"""
    def test_unit_interval(self):
        """lists of random values in [0, 1] chosen. Minimum distance cannot be greater than 1."""
        n = numpy.random.rand(5, 3)
        m = numpy.random.rand(10, 3)
        min_dist = isambard.geometry.closest_distance(points1=n, points2=m)
        self.assertLessEqual(min_dist, 1.0)

    def test_zero_case(self):
        """ Random arrays with at least one shared element should have a closest distance of 0."""
        n = numpy.random.rand(5, 3)
        m = numpy.random.rand(10, 3)
        # force first elements of otherwise random arrays to be equal.
        m[0] = n[0]
        min_dist = isambard.geometry.closest_distance(points1=n, points2=m)
        self.assertAlmostEqual(min_dist, 0)


class MinimalDistanceBetweenLines(unittest.TestCase):
    """ tests for minimal_distance_between_lines. """
    def setUp(self):
        self.vectors = random_vectors(n=100, vector_length=3)

    def test_simple_example(self):
        """ simple example taken from http://www.math.jhu.edu/~js/Math202/example1.pdf """
        a = numpy.array([-2, 1, -1])
        b = a + numpy.array([2, 3, -1])
        c = numpy.array([1, -1, 2])
        d = c + numpy.array([-1, 2, 4])
        self.assertTrue(numpy.isclose(isambard.geometry.minimal_distance_between_lines(a, b, c, d, segments=False),
                                      11.0/numpy.sqrt(6.0)))

    def test_minimal_line_distance(self):
        for i in range(100):
            a, b, c, d = random.sample(self.vectors, 4)
            segments = random.choice([True, False])
            d1 = isambard.geometry.minimal_distance_between_lines(a, b, c, d, segments=segments)
            d2 = isambard.geometry.distance(a, c)
            d3 = isambard.geometry.distance(a, d)
            d4 = isambard.geometry.distance(b, c)
            d5 = isambard.geometry.distance(b, d)
            at = numpy.isclose(min([d1, d2, d3, d4, d5]), d1) or (numpy.isclose(d1, 0, atol=1e-07))
            self.assertTrue(at)




class IntersectionStartAndEndTestCase(unittest.TestCase):
    def test_example_in_function_notes(self):
        a = [numpy.array([0, 0, x]) for x in range(11)]
        b = [numpy.array([10, 0, x]) for x in numpy.linspace(7, 12, num=1000)]
        self.assertEqual(isambard.geometry.intersection_start_and_end(subject=a, reference=b), (7, 10))
        self.assertEqual(isambard.geometry.intersection_start_and_end(subject=b, reference=a), (0, 599))


class FindTransformationsTestCase(unittest.TestCase):
    """Tests for isambard.tools.geometry.find_transformations"""
    def setUp(self):
        self.vectors = [numpy.array([random.randint(-100, 100) * random.random()
                                     for _ in range(3)]) for _ in range(100)]

    def test_translation(self):
        # translating s1 should give s2.
        at = []
        for _ in range(100):
            a, b, c, d = random.sample(self.vectors, 4)
            translation = isambard.geometry.find_transformations(s1=a, e1=b, s2=c, e2=d)[0]
            at.append(numpy.allclose(a + translation, c))
        self.assertTrue(all(at))

    def test_rotation(self):
        at = []
        for _ in range(100):
            a, b, c, d = random.sample(self.vectors, 4)
            angle, axis, point = isambard.geometry.find_transformations(a, b, c, d)[1:]
            q = isambard.geometry.Quaternion.angle_and_axis(angle, axis)
            t = q.rotate_vector(b, point=point)
            x = isambard.geometry.angle_between_vectors(t - a, d - c)
            at.append(numpy.isclose(x, 0.0, atol=0.0001))
        self.assertTrue(all(at))

    def test_axis_not_origin(self):
        at = []
        for _ in range(100):
            a, b, c, d = random.sample(self.vectors, 4)
            axis = isambard.geometry.find_transformations(a, b, c, d)[2]
            at.append(not numpy.allclose(axis, [0, 0, 0]))
        self.assertTrue(all(at))


class CoordinateSystemTestCase(unittest.TestCase):
    """Tests for coordinate system conversion functions in isambard.tools.geometry"""
    def setUp(self):
        n = 100
        # non-zero angles for azimuths and zeniths.
        self.azimuths = []
        azimuths = random_angles(n=n)
        for a in azimuths:
            if a == 0:
                self.azimuths.append(1)
            else:
                self.azimuths.append(a)
        self.zeniths = []
        zeniths = random_angles(n=n)
        for z in zeniths:
            if z == 0:
                self.zeniths.append(1)
            else:
                self.zeniths.append(z)
        # positive values for radii
        self.radii = random_floats(n=n, min_val=1, max_val=100)
        self.cartesians = random_vectors(n=n)

    def test_spherical_to_cartesian_and_back(self):
        """ Convert spherical coordinates to cartesian and back again,
         check returned values are the same as initial values."""
        to_cartesian = [isambard.geometry.spherical_to_cartesian(r, a, z) for r, a, z in zip(
            self.radii, self.azimuths, self.zeniths)]
        back_to_spherical = [isambard.geometry.cartesian_to_spherical(x, y, z) for x, y, z in to_cartesian]
        a1 = numpy.allclose(self.radii, [x[0] for x in back_to_spherical])
        a2 = numpy.allclose(self.azimuths, [x[1] for x in back_to_spherical])
        a3 = numpy.allclose(self.zeniths, [x[2] for x in back_to_spherical])
        at = all([a1, a2, a3])
        self.assertTrue(at)

    def test_cartesian_to_spherical_and_back(self):
        """ Convert cartesian coordinates to spherical and back again,
         check returned values are the same as initial values."""
        to_spherical = [isambard.geometry.cartesian_to_spherical(x, y, z) for x, y, z in self.cartesians]
        back_to_cartesian = [isambard.geometry.spherical_to_cartesian(r, a, z) for r, a, z in to_spherical]
        at = numpy.allclose(self.cartesians, back_to_cartesian)
        self.assertTrue(at)

    def test_cylindrical_to_cartesian_and_back(self):
        """ Convert cylindrical coordinates to cartesian and back again,
         check returned values are the same as initial values."""
        to_cartesian = [isambard.geometry.cylindrical_to_cartesian(r, a, z)
                        for r, a, z in zip(self.radii, self.azimuths, [v[2] for v in self.cartesians])]
        back_to_cylindrical = [isambard.geometry.cartesian_to_cylindrical(x, y, z) for x, y, z in to_cartesian]
        a1 = numpy.allclose(self.radii, [x[0] for x in back_to_cylindrical])
        a2 = numpy.allclose(self.azimuths, [x[1] for x in back_to_cylindrical])
        a3 = numpy.allclose([v[2] for v in self.cartesians], [x[2] for x in back_to_cylindrical])
        at = all([a1, a2, a3])
        self.assertTrue(at)

    def test_cartesian_to_cylindrical_and_back(self):
        """ Convert cartesian coordinates to cylindrical and back again,
         check returned values are the same as initial values."""
        to_cylindrical = [isambard.geometry.cartesian_to_cylindrical(x, y, z) for x, y, z in self.cartesians]
        back_to_cartesian = [isambard.geometry.cylindrical_to_cartesian(r, a, z) for r, a, z in to_cylindrical]
        at = numpy.allclose(self.cartesians, back_to_cartesian)
        self.assertTrue(at)


class QuaternionTestCase(unittest.TestCase):
    """Tests for isambard.tools.quaternion.Quaternion class"""
    def setUp(self):
        """ Generate 100 random Quaternions for testing"""
        self.quats = []
        for x in range(100):
            r, i, j, k = numpy.random.rand(4) * [numpy.random.choice(range(-100, 100)) for _ in range(4)]
            q = isambard.geometry.Quaternion(r=r, i=i, j=j, k=k)
            self.quats.append(q)

        # Random angles in range [0, 180]
        angles = [numpy.rad2deg(numpy.arccos(numpy.cos(q.r))) for q in self.quats]
        angles = [x for x in angles if (x > 0) and (x < 180)]
        self.angles = angles

        self.vectors = [numpy.random.rand(3) for _ in range(len(self.angles))]
        return

    def test_norm_of_rotation_is_1(self):
        """ q.as_rotation() should always be a unit quaternion """
        desired = numpy.array([1.0] * len(self.quats))
        actual = [q.as_rotation().norm for q in self.quats]
        at = numpy.allclose(actual, desired)
        self.assertTrue(at)

    def test_rotation_conjugate_is_inverse(self):
        """ Unit quaternions should satisfy q.conjugate == ~q """
        conjugates = [q.as_rotation().conjugate.as_array for q in self.quats]
        inverses = [(~q.as_rotation()).as_array for q in self.quats]
        at = numpy.allclose(conjugates, inverses)
        self.assertTrue(at)

    def test_norm_and_division(self):
        """ Dividing a quaternion by its norm should always result in a unit quaternion """
        quats = [q for q in self.quats if q.norm != 0]
        desired = numpy.array([1.0] * len(quats))
        actual = [(q / q.norm).norm for q in quats]
        at = numpy.allclose(actual, desired)
        self.assertTrue(at)

    def test_both_multiplications_with_floats(self):
        """ Left- and right- multiplications by floats should be identical """
        results = []
        for q in self.quats:
            x = random.random()
            results.append(numpy.allclose((x * q).as_array, (q * x).as_array))
        at = all(results)
        self.assertTrue(at)

    def test_rotation_by_zero_is_unit(self):
        """ Applying conjugation to any unit quaternion with zero real part should return Q([1, 0, 0, 0]) """
        angle = 0
        desired = numpy.array([1, 0, 0, 0])
        vectors = self.vectors
        q_vectors = [isambard.geometry.Quaternion.real_and_vector(angle, v).as_rotation() for v in vectors]
        actual = [q.apply_conjugation(vq).as_array for q, vq in zip(self.quats, q_vectors)]
        at = numpy.allclose(actual, desired)
        self.assertTrue(at)

    def test_rotation_and_extraction(self):
        """ Using q.as_rotation() followed by q.extract_angle_and_axis() should get you back to the start """
        angles = self.angles
        # Random taken from self.quats
        vectors = [q.vector for q in self.quats[:len(angles)]]
        qs = [isambard.geometry.Quaternion.real_and_vector(angle, v).as_rotation() for angle, v in zip(angles, vectors)]
        r_angles, r_vectors = zip(*[q.extract_angle_and_axis(radians=False) for q in qs])
        # returned vectors should be unit_vectors of initial vectors.
        vectors = [isambard.geometry.unit_vector(v) for v in vectors]
        at1 = numpy.allclose(angles, r_angles)
        at2 = numpy.allclose(vectors, r_vectors)
        at = at1 and at2
        self.assertTrue(at)

    def test_angle_and_axis(self):
        desired = [isambard.geometry.Quaternion.real_and_vector(a, v).as_rotation().as_array for a, v in zip(
            self.angles, self.vectors)]
        actual = [isambard.geometry.Quaternion.angle_and_axis(a, v).as_array for a, v in zip(self.angles, self.vectors)]
        at = numpy.allclose(desired, actual)
        self.assertTrue(at)

    def test_apply_conjugation_inputs(self):
        """ apply_conjugation should work with vectors or quaternions as input
        """
        rot_quats = [isambard.geometry.Quaternion.angle_and_axis(a, v) for a, v in zip(self.angles, self.vectors)]
        test_vectors = [numpy.random.rand(3) for _ in range(len(self.angles))]
        test_quats = [isambard.geometry.Quaternion.real_and_vector(0, v) for v in test_vectors]
        q_results = [r.apply_conjugation(x).vector for r, x in zip(rot_quats, test_quats)]
        v_results = [r.apply_conjugation(x) for r, x in zip(rot_quats, test_vectors)]
        at = numpy.allclose(q_results, v_results)
        self.assertTrue(at)


class AxisTestCase(unittest.TestCase):
    """Tests for isambard.tools.geometry.Axis class"""
    def setUp(self):
        """ Generate n random starts and ends for testing"""
        n = 100
        starts = random_vectors(n=n)
        ends = random_vectors(n=n)
        self.axes = [isambard.geometry.Axis(start=s, end=e) for s, e in zip(starts, ends)]

    def test_units_perpendicular(self):
        desired = [90.0] * len(self.axes)
        unit_to_rad = [isambard.geometry.angle_between_vectors(a.unit_tangent, a.unit_normal) for a in self.axes]
        unit_to_tan = [isambard.geometry.angle_between_vectors(a.unit_tangent, a.unit_binormal) for a in self.axes]
        rad_to_tan = [isambard.geometry.angle_between_vectors(a.unit_normal, a.unit_binormal) for a in self.axes]
        a1 = numpy.allclose(unit_to_rad, desired)
        a2 = numpy.allclose(unit_to_tan, desired)
        a3 = numpy.allclose(rad_to_tan, desired)
        at = a1 and a2 and a3
        self.assertTrue(at)

    def test_ax_unit_and_length(self):
        end_test = [a.start + (a.length * a.unit_tangent) for a in self.axes]
        actual_ends = [a.end for a in self.axes]
        at = numpy.allclose(end_test, actual_ends)
        self.assertTrue(at)


class HelicalCurveTestCase(unittest.TestCase):
    """Tests for isambard.ampal.specifications.primitive_paths.HelicalCurve"""
    def setUp(self):
        n = 100
        self.alphas = random_angles(n=n, min_val=0, max_val=90)
        self.radii = [random.random() * 100 for _ in range(n)]
        self.pitches = [random.random() * 1000 for _ in range(n)]
        self.curves = [isambard.geometry.HelicalCurve.alpha_and_pitch(
            alpha=a, pitch=p) for a, p in zip(self.alphas, self.pitches)]
        return

    def test_instantiation_pitch_and_radius(self):
        helical_curves = [isambard.geometry.HelicalCurve.pitch_and_radius(
            pitch=p, radius=r) for p, r in zip(self.pitches, self.radii)]
        a1 = numpy.allclose([h.pitch for h in helical_curves], self.pitches)
        a2 = numpy.allclose([h.radius for h in helical_curves], self.radii)
        self.assertTrue(a1 and a2)

    def test_instantiation_alpha_and_radius(self):
        helical_curves = [isambard.geometry.HelicalCurve.alpha_and_radius(
            alpha=a, radius=r) for a, r in zip(self.alphas, self.radii)]
        a1 = numpy.allclose([h.alpha for h in helical_curves], self.alphas)
        a2 = numpy.allclose([h.radius for h in helical_curves], self.radii)
        self.assertTrue(a1 and a2)

    def test_instantiation_alpha_and_pitch(self):
        helical_curves = [isambard.geometry.HelicalCurve.alpha_and_pitch(
            alpha=a, pitch=p) for a, p in zip(self.alphas, self.pitches)]
        a1 = numpy.allclose([h.alpha for h in helical_curves], self.alphas)
        a2 = numpy.allclose([h.pitch for h in helical_curves], self.pitches)
        self.assertTrue(a1 and a2)

    def test_t_from_arc_length(self):
        """ t_from_arc_length and arc_length should be inverses of each other. """
        t_values = [random.random() * random.choice(range(-1000, 1000)) for _ in range(len(self.alphas))]
        calculated_t_values = [hc.t_from_arc_length(hc.arc_length(x)) for hc, x in zip(self.curves, t_values)]
        at = numpy.allclose(t_values, calculated_t_values)
        self.assertTrue(at)

    def test_get_coords_n_points(self):
        n_point_vals = [random.randint(1, 100) for _ in range(100)]
        len_coords = [len(hc.get_coords(n_points=x)) for hc, x in zip(self.curves, n_point_vals)]
        at = numpy.allclose(n_point_vals, len_coords)
        self.assertTrue(at)


__author__ = 'Christopher W. Wood, Jack W. Heal'


