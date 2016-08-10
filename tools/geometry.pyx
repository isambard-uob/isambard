"""Cythonised versions of basic geometric operations that ISAMBARD requires.

If Cython is not detected the geometry module falls back to the Python
implementation. As a result, all functions that are available to the user MUST
return the same type the Python version of the code.

Functions that have names preceded by an underscore indicate that they are not
to be used by users. However, all of these functions will generally be defined
using cdef and so will not be imported anyway.
"""

#cython: embedsignature=True

from libc.math cimport acos, cos, sin, atan2, sqrt, pow, fmax, fmin, M_PI
from numbers import Number

import numpy
import itertools
import sys

cdef double SMALL_NUMBER = 1e-7
cdef double SMALLEST_FLOAT = sys.float_info.min

numpy.seterr(divide='ignore', invalid='ignore')  # Suppress invalid warnings, as these are caught anyway.

# *******************
# Cythonised geometry
# *******************



cdef double _deg_to_rad(double angle):
    return (angle/180.0) * M_PI


cdef double _rad_to_deg(double angle):
    return (angle/M_PI) * 180.0


cdef double _dot_product_2d(double[2] a, double[2] b):
    cdef double dot
    dot = (a[0] * b[0]) + (a[1] * b[1])
    return dot


cdef double _dot_product_3d(double[3] a, double[3] b):
    cdef double dot
    dot = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])
    return dot


cdef _cross_product(double[3] a, double[3] b):
    cdef double cross[3]
    cross[0] = a[1]*b[2] - a[2]*b[1]
    cross[1] = a[2]*b[0] - a[0]*b[2]
    cross[2] = a[0]*b[1] - a[1]*b[0]
    return cross


cdef _subtract_vectors_3d(double[3] a, double[3] b):
    cdef double subtract[3]
    subtract[0] = a[0] - b[0]
    subtract[1] = a[1] - b[1]
    subtract[2] = a[2] - b[2]
    return subtract


cdef double _magnitude(double[3] a):
    cdef double mag
    mag = sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2))
    return mag


cpdef double distance(a, b):
    """Calculates the Euclidian distance between two points in R^3.

    Parameters
    ----------
    a : list or tuple or numpy.array
    b : list or tuple or numpy.array

    Returns
    -------
    dist : double
        Euclidean distance between a and b.
    """
    cdef double a_x, a_y, a_z, b_x, b_y, b_z, dist
    a_x, a_y, a_z = a
    b_x, b_y, b_z = b
    dist = sqrt(pow((b_x - a_x), 2) + pow((b_y - a_y), 2) + pow((b_z - a_z), 2))
    return dist


def unit_vector(a):
    """Return vector normalised by length.

    Parameters
    ----------
    a : list or tuple or numpy.array
        The vector to be normalised.

    Returns
    -------
    v_out : numpy.array
        A numpy.array that is a unit vector.

    Raises
    ------
    ZeroDivisionError
        If the vector has length 0.
    """
    cdef double a_x, a_y, a_z, vector_length
    cdef double[3] a_in = a
    vector_length = _magnitude(a_in)
    if vector_length != 0:
        v_out = [a_in[0]/vector_length, a_in[1]/vector_length, a_in[2]/vector_length]
        return numpy.array(v_out)
    else:
        raise ZeroDivisionError("Vector must be of non-zero length.")


cdef _unit_vector(double[3] a):
    cdef double[3] v_out
    vector_length = _magnitude(a)
    v_out[0] = a[0]/vector_length
    v_out[1] = a[1]/vector_length
    v_out[2] = a[2]/vector_length
    return v_out


cpdef double angle_between_vectors(a, b, radians=False):
    """Calculate angle between two vectors using the dot product.

    Notes
    -----
    All angle-based functions in tools_geometry will return values in degrees by default.
    The radians flag can be used to return values in radians.
    If a or b is very close to the origin (norm < 1e-7), returns 0.0.

    Parameters
    ----------
    a : list or tuple or numpy.array.
    b : list or tuple or numpy.array.
    radians : bool, optional
        If True, returns angle in radians.
        If False, returns angle in degrees.
        Default value is False.

    Returns
    -------
    angle : double
        The angle between vectors v1 and v2.
    """
    cdef double[3] a_v, b_v
    cdef double mag_a, mag_b, angle, cos_angle
    a_v = a
    b_v = b
    mag_a = _magnitude(a_v)
    mag_b = _magnitude(b_v)
    if (mag_a < SMALL_NUMBER) or (mag_b < SMALL_NUMBER):
        return 0.0
    cos_angle = _dot_product_3d(a_v, b_v) / (mag_a * mag_b)
    # take care of round-off error.
    cos_angle = fmin(cos_angle, 1.0)
    cos_angle = fmax(-1.0, cos_angle)
    angle = acos(cos_angle)
    if radians:
        return angle
    else:
        return _rad_to_deg(angle)

cpdef int is_acute(a, b):
    """ Return True if angle between two vectors is less than 90 degrees.

    Parameters
    ----------
    a : list or tuple or numpy.array
    b : list or tuple or numpy.array

    Returns
    -------
    bool
        True if angle between v1 and v2 is less than 90 degrees.
        False otherwise.

    """
    cdef double angle = angle_between_vectors(a, b)
    return 1 if angle < 90.0 else 0


cpdef double dihedral(a, b, c, d, radians=False):
    """Calculates the dihedral angle defined by four points in R^3.

    The points a, b and c form a plane in R^3.
    The points b, c and d form a plane in R^3.
    The dihedral angle is the angle between these two planes.

    Notes
    -----
    See http://en.wikipedia.org/wiki/Dihedral_angle and https://en.wikipedia.org/wiki/Atan2 for more information.
    Dihedral between 4 points.
    a, b, c, d are all 3-vectors (x,y,z) coordinates.
    (a,b,c) and (b,c,d) must both be non-colinear.

    Parameters
    ----------
    a : list or tuple or numpy.array
    b : list or tuple or numpy.array
    c : list or tuple or numpy.array
    d : list or tuple or numpy.array
    radians : bool, optional
        If True, returns angle in radians.
        If False, returns angle in degrees.
        Default value is False.

    Returns
    -------
    angle : double
        The dihedral angle formed between the points a, b, c and d.
    """
    cdef double[3] a_v, b_v, c_v, d_v, p, q, r, n1, n2, n3
    cdef double dist
    a_v = a
    b_v = b
    c_v = c
    d_v = d
    # cross product between b and d is perpendicular to both v1 and v2.
    p = _subtract_vectors_3d(b_v, a_v)
    q = _subtract_vectors_3d(c_v, b_v)
    r = _subtract_vectors_3d(d_v, c_v)
    n1 = _cross_product(p, q)
    n2 = _cross_product(q, r)
    if (_magnitude(n1) < SMALL_NUMBER) or (_magnitude(n2) < SMALL_NUMBER):
        return 0.0
    dihedral_angle = angle_between_vectors(n1, n2, radians=radians)
    n3 = _cross_product(n1, n2)
    if (not _magnitude(q) < SMALL_NUMBER) and (not _magnitude(n3) < SMALL_NUMBER):
        if angle_between_vectors(q, n3, radians=True) > SMALL_NUMBER:
            dihedral_angle = -dihedral_angle
    return dihedral_angle


def find_foot(a, b, p):
    """Finds the perpendicular foot of a point (p) and a line (ab).

    The points a and b define a line in R^3.
    There is a point on this line, q, that has the shortest Euclidian distance from the point p.
    This point is the perpendicular foot, and is returned by this function.

    Notes
    -----
    See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html for a useful work-through of the maths.
    Parameterise the line ab by t, so that a point on the line is of the form: a + (b-a)*t.
    Define f(t) = [distance(t, p)]^2. Want t such that f(t) is at a minimum.
    Solve f'(t) = 0 for t.
    Substituting this value of t into a + (b-a)*t yields the desired point q.

    Parameters
    ----------
    a : list or tuple or numpy.array
    b : list or tuple or numpy.array
        Together, the points a and b define a line in R^3.
    p : list or tuple or numpy.array
        A point in R^3.

    Returns
    -------
    q : list or tuple or numpy.array
        The point on the line between a and b that is closest to q.
    """
    cdef double a_x, a_y, a_z, b_x, b_y, b_z, p_x, p_y, p_z
    cdef double[3] c, d, q

    a_x, a_y, a_z = a
    b_x, b_y, b_z = b
    p_x, p_y, p_z = p

    c[0] = b_x - a_x
    c[1] = b_y - a_y
    c[2] = b_z - a_z

    d[0] = a_x - p_x
    d[1] = a_y - p_y
    d[2] = a_z - p_z

    t = -(_dot_product_3d(d, c)) / (pow(b_x - a_x, 2) + pow(b_y - a_y, 2) + pow(b_z - a_z, 2))
    q[0] = a_x + c[0] * t
    q[1] = a_y + c[1] * t
    q[2] = a_z + c[2] * t
    return numpy.array(q)


def find_foot_on_plane(a, b, c, p):
    """Finds the perpendicular foot of a point (p) and a plane defined by points a, b and c.

    The points a, b and c must define a plane in R^3 (i.e. they cannot be collinear).
    There is a point on this plane, q, that has the shortest Euclidian distance from the point p.
    This point is the perpendicular foot, and is returned by this function.

    Notes
    -----
    This page:
     http://math.stackexchange.com/questions/723937/find-the-point-on-a-plane-3x-4y-z-1-that-is-closest-to-1-0-1
     gives a work-through of this kind of problem for a specific example.
    The solution uses the Lagrangian multiplier: https://en.wikipedia.org/wiki/Lagrange_multiplier
    We look for the minimum values of the distance function f (squared to make it easier!)
     f(x, y, z) = (x - p[0])^2 + (y - p[1])^2 + (z - p[2])^2
    subject to the constraint g(x, y, z) = 0, where g is the equation of the plane.
    The Langrangian multiplier L = f(x, y, z) - lambda*(g(x, y, z)).
    We calculate the partial deriviatves of this (wrt x, y, z and lambda).
    Solve first for lambda gives the equation for the lagrangian_lambda value in this function.
    Then substitute lambda back in to get the solution for (x, y, z) = q.

    Parameters
    ----------
    a : list or tuple or numpy.array
    b : list or tuple or numpy.array
    c : list or tuple or numpy.array
        Together, the points a, b and c define a plane in R^3.
    p : list or tuple or numpy.array
        A point in R^3.

    Returns
    -------
    s : numpy.array
        The point on the plane that is closest to q.
    """
    cdef double[3] aa, ba, ca, pa, q, r, s
    aa = a
    ba = b
    ca = c
    pa = p

    q[0] = ba[0] - aa[0]
    q[1] = ba[1] - aa[1]
    q[2] = ba[2] - aa[2]

    r[0] = ca[0] - ba[0]
    r[1] = ca[1] - ba[1]
    r[2] = ca[2] - ba[2]
    # n is normal to the plane defined by a, b and c.
    cdef double[3] n = _cross_product(q, r)

    if _magnitude(n) == 0:
        raise ValueError('The points a, b and c are collinear and therefore do not define a plane.')

    cdef double d = _dot_product_3d(n, aa)
    # plane equation is numpy.dot(n, [x, y, z]) - d = 0.

    half_lagrangian_lambda = (_dot_product_3d(n, pa) - d) / (_dot_product_3d(n, n))

    s[0] = p[0] - (half_lagrangian_lambda * n[0])
    s[1] = p[1] - (half_lagrangian_lambda * n[1])
    s[2] = p[2] - (half_lagrangian_lambda * n[2])

    return numpy.array(s)


def radius_of_circumcircle(a, b, c):
    """Returns the radius the circumcircle of a triangle, when supplied with three points (p, q, r).

    Notes
    -----
    See http://www.analyzemath.com/Geometry/Circumcircle/CircumcircleRadius.html
    for more information.
    Three (non-colinear) points define a traingle.
    The circumcircle is the circle that passes though these points.

    Parameters
    ----------
    p : list or tuple or numpy.array
    q : list or tuple or numpy.array
    r : list or tuple or numpy.array

    Returns
    -------
    radius : float
    """
    # calculate edge lengths of the triangle.
    cdef double a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z, p, q, r, cos_angle, angle, radius

    a_x, a_y, a_z = a
    b_x, b_y, b_z = b
    c_x, c_y, c_z = c

    p = pow((b_x - c_x), 2) + pow((b_y - c_y), 2) + pow((b_z - c_z), 2)
    q = pow((a_x - c_x), 2) + pow((a_y - c_y), 2) + pow((a_z - c_z), 2)
    r = pow((a_x - b_x), 2) + pow((a_y - b_y), 2) + pow((a_z - b_z), 2)

    # cosine rule
    cos_angle = (q + r - p) / (2 * sqrt(q) * sqrt(r))

    # Use side length and angle to calculate the circumcircle radius.
    angle = acos(cos_angle)
    radius = sqrt(p) / (2 * sin(angle))

    return radius

# TODO Remove this function and its tests when find_and_fit_loops no longer uses it.
def rotation_matrix(angle, unit_vect, point=None):
    """Find rotation matrix for rotating about unit_vect by angle.

    Notes
    -----
    See https://en.wikipedia.org/wiki/Rotation_matrix for more information.
    In particular, the section headed 'Rotation matrix from axis and angle'.

    Parameters
    ----------
    angle : float
        Angle about which to rotate.
    unit_vect : list or tuple or numpy.array
        Vector about which to rotate.
    point : list or tuple or numpy.array, optional.
        Specify a point about which to rotate.
        If point is None, then the rotation will be about the origin.
        Default value is None.

    Returns
    -------
    m : numpy.array
        A (4 x 4) rotation matrix.
    """
    cdef double sina, cosa
    # TODO Confirm whether we want the input angle to be in degrees or radians. Change unit vector to internally do this.
    sina = sin(angle)
    cosa = cos(angle)
    # TODO use the unit_vector function here so that the function can take a vector of arbitrary length.
    direction = numpy.array([float(x) for x in unit_vect])

    # rotation matrix around unit vector
    r = numpy.diag([cosa, cosa, cosa])
    r += numpy.outer(direction, direction) * (1.0 - cosa)

    # cross product matrix of the unit_vector multiplied by sina
    direction *= sina
    r += numpy.array([[0.0, -direction[2], direction[1]],
                      [direction[2], 0.0, -direction[0]],
                      [-direction[1], direction[0], 0.0]])
    m = numpy.identity(4)
    m[:3, :3] = r
    if point is not None:
        # rotation not around origin
        point = numpy.array(point[:3], dtype=numpy.float64, copy=False)
        m[:3, 3] = point - numpy.dot(r, point)
    return m


def apply_unit_quaternion(vector, u_in, double s):
    """Applies the rotation defined by the unit quaternion to a 3D vector.

    Parameters
    ----------
    vector : [float, float, float]
        3D coordinate.
    u_in : double[3]
        Imaginary component of the unit quaternion.
    s : double
        Real component of the quaternion.

    Returns
    -------
    r : double[3]
        3D coordinate.
    """
    cdef double[3] r, u, v, c_uv
    u = u_in
    v = vector
    cdef double d_uv = _dot_product_3d(u, v) * 2
    cdef double d_uu = _dot_product_3d(u, u)
    c_uv = _cross_product(u, v)
    cdef double stp2 = pow(s, 2) - d_uu
    cdef double s2 = s * 2
    r[0] = (d_uv * u[0]) + (stp2 * v[0]) + (s2 * c_uv[0])
    r[1] = (d_uv * u[1]) + (stp2 * v[1]) + (s2 * c_uv[1])
    r[2] = (d_uv * u[2]) + (stp2 * v[2]) + (s2 * c_uv[2])
    return numpy.array(r)


def minimal_distance_between_lines(p0, p1, q0, q1, segments=True):
    """Finds the minimal distance between two the two lines joining p0 to p1 and q0 to q1.

    The points p0, p1, q0 and d are used to construct two lines in R^3.
    The first line L1 is from p0 to p1, the second line L2 is from q0 to q1.
    The minimal distance between these lines is the length of the line that
    connects them and is perpendicular to both L1 and L2.

    Notes
    -----
    See https://en.wikipedia.org/wiki/Skew_lines#Distance for more information.
    Also http://geomalgorithms.com/a07-_distance.html

    Parameters
    ----------
    p0 : list or tuple or numpy.array
    p1 : list or tuple or numpy.array
    q0 : list or tuple or numpy.array
    q1 : list or tuple or numpy.array
    segments : bool
        If True, the points returned must lie on line segments p0->p1 and q0->q1, rather than their infinite extensions.

    Returns
    -------
    dist : float
        The distance between the two lines.
    """
    return distance(*closest_points_on_lines(p0, p1, q0, q1, segments=segments))


def closest_points_on_lines(p0, p1, q0, q1, segments=True):
    """
    Notes
    -----
    p0, p1, q0, q1 describe two lines, L1 (p0 -> p1) and L2 (q0 -> q1).
    L1 (2) is the infinite extension of the line that passes through a (q0) and has direction p0->p1 (q0->q1).
    Two points point_1, and point_2 on L1 and L2 respectively, are returned.
    The distance between point_1 and point_2 is smaller than any other points on L1 and L2.
    Code adapted from http://geomalgorithms.com/a07-_distance.html, where a more detailed explanation is available.


    Parameters
    ----------
    p0 : list or tuple or numpy.array
    p1 : list or tuple or numpy.array
    q0 : list or tuple or numpy.array
    q1 : list or tuple or numpy.array
    segments : bool
        If True, the points returned must lie on line segments p0->p1 and q0->q1, rather than their infinite extensions.

    Returns
    -------
    point_1: numpy.array
    point_2: numpy.array

    """
    cdef double[3] p0_v, p1_v, q0_v, q1_v, u, v, w0, point_1, point_2
    cdef double a, b, c, d, e, denom, tc, tn, td, sc, sn, sd
    p0_v = p0
    p1_v = p1
    q0_v = q0
    q1_v = q1
    u = _subtract_vectors_3d(p1_v, p0_v)
    v = _subtract_vectors_3d(q1_v, q0_v)
    w0 = _subtract_vectors_3d(p0_v, q0_v)
    a = _dot_product_3d(u, u)
    b = _dot_product_3d(u, v)
    c = _dot_product_3d(v, v)
    d = _dot_product_3d(u, w0)
    e = _dot_product_3d(v, w0)
    denom = a * c - pow(b, 2)
    # if both lines close to zero in length, return distance closest end points.
    if a < SMALL_NUMBER and c < SMALL_NUMBER:
        np1 = p0
        np2 = q0
        for x1, x2 in itertools.product([p0, p1], [q0, q1]):
            d = distance(x1, x2)
            if d < distance(np1, np2):
                np1 = x1
                np2 = x2
        return np2, np1
    # if line1 close to zero in length (but not line2), return any point on line1 and its foot on line2
    elif a < SMALL_NUMBER:
        p = find_foot(q0, q1, p0)
        return p0, p
    # if line2 close to zero in length (but not line1), return any point on line2 and its foot on line1
    elif c < SMALL_NUMBER:
        p = find_foot(p0, p1, q0)
        return q0, p
    if not segments:
        if denom < SMALL_NUMBER:
            sc = 0.0
            if b > c:
                tc = d / b
            else:
                tc = e / c
        else:
            sc = (b*e - c*d) / denom
            tc = (a*e - b*d) / denom
    else:
        sc = denom
        sn = denom
        sd = denom
        tc = denom
        tn = denom
        td = denom
        if denom < SMALL_NUMBER:
            sn = 0.0
            sd = 1.0
            tn = e
            td = c
        else:
            sn = ((b * e) - (c * d))
            tn = ((a * e) - (b * d))
            if sn < 0.0:
                sn = 0.0
                tn = e
                td = c
            elif sn > sd:
                sn = sd
                tn = e + b
                td = c
        if tn < 0.0:
            tn = 0.0
            if -d < 0.0:
                sn = 0.0
            elif -d > a:
                sn = sd
            else:
                sn = -d
                sd = a
        elif tn > td:
            tn = td
            if (-d + b) < 0.0:
                sn = 0.0
            elif (-d + b) > a:
                sn = sd
            else:
                sn = (-d + b)
                sd = a
        sc = 0 if sn < SMALL_NUMBER else sn / sd
        tc = 0 if tn < SMALL_NUMBER else tn / td
    point_1[0] = p0_v[0] + (sc * u[0])
    point_1[1] = p0_v[1] + (sc * u[1])
    point_1[2] = p0_v[2] + (sc * u[2])
    point_2[0] = q0_v[0] + (tc * v[0])
    point_2[1] = q0_v[1] + (tc * v[1])
    point_2[2] = q0_v[2] + (tc * v[2])
    return point_1, point_2



# *********
# Coordinate system conversion functions.
# *********
def spherical_to_cartesian(double radius, double azimuth, double zenith, radians=False):
    """ Converts from spherical to cartesian coordinates system

    Notes
    -----
    Spherical coordinates are here defined by a radius, an azimuth angle and a zenith angle.
    The zenith is the angle between the point and the z-axis.
    So: zenith = 0 => lies on z-axis, zenith = 90 => lies in x-y plane, zenith = 180 => lies antiparallel to z-axis.
    The azimuth is the angle in the x-y plane, with counter-clockwise direction being positive.
    So: azimuth = 0 => lies on x-axis, azimuth = 90 => lies in y axis, zenith = 180 => lies antiparallel to x-axis.
    See https://en.wikipedia.org/wiki/Spherical_coordinate_system
    and/or http://mathworld.wolfram.com/SphericalCoordinates.html
    for more information.
    If zenith = 0, then x and y are both 0 regardless of azimuth value.
    If azimuth = 0, then y is 0 regardless of zenith value.

    Parameters
    ----------
    radius : float
    azimuth : float
    zenith : float
    radians : bool
        True if azimuth and zenith are given in radians instead of degrees.

    Returns
    -------
    numpy.array([x, y, z])
        cartesian coordinates corresponding to the input spherical coordinates.
    """
    if not radians:
        azimuth = _deg_to_rad(azimuth)
        zenith = _deg_to_rad(zenith)
    cdef double x, y, z
    x = radius * sin(zenith) * cos(azimuth)
    y = radius * sin(zenith) * sin(azimuth)
    z = radius * cos(zenith)
    return numpy.array([x, y, z])

def cartesian_to_spherical(double x, double y, double z, radians=False):
    """ Converts from cartesian to spherical coordinates system

    Notes
    -----
    See spherical_to_cartesian for more information.

    Parameters
    ----------
    x : float
    y : float
    z : float
    radians : bool
        True if returned azimuth and zenith are in radians instead of degrees.

    Returns
    -------
    radius : float
    azimuth : float
    zenith : float
        Spherical coordinate parameters. See spherical_to_cartesian for more information.
    """
    cdef double radius, azimuth, zenith
    if (x == 0) and (y == 0) and (z == 0):
        return 0.0, 0.0, 0.0
    radius = sqrt(x**2 + y**2 + z**2)
    azimuth = atan2(y, x)
    zenith = acos(z / radius)
    if not radians:
        azimuth = _rad_to_deg(azimuth)
        zenith = _rad_to_deg(zenith)
    return radius, azimuth, zenith

def cylindrical_to_cartesian(double radius, double azimuth, double z, radians=False):
    """ Converts from cylindrical to cartesian coordinates system

    Notes
    -----
    Cylindrical coordinates are here defined by a radius, an azimuth angle and a z value.
    The azimuth is the angle in the x-y plane, with counter-clockwise direction being positive.
    So: azimuth = 0 => lies on x-axis, azimuth = 90 => lies in y axis, zenith = 180 => lies antiparallel to x-axis.
    The radius and the z value are the radius and height of the cylinder, respectively.
    See https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
    and/or http://mathworld.wolfram.com/CylindricalCoordinates.html
    for more information.

    Parameters
    ----------
    radius : float
    azimuth : float
    z : float
    radians : bool
        True if azimuth is given in radians instead of degrees.

    Returns
    -------
    numpy.array([x, y, z])
        cartesian coordinates corresponding to the input cylindrical coordinates.
    """
    cdef double x, y
    if not radians:
        azimuth = _deg_to_rad(azimuth)
    x = radius * cos(azimuth)
    y = radius * sin(azimuth)
    return numpy.array([x, y, z])

def cartesian_to_cylindrical(double x, double y, double z, radians=False):
    """ Converts from cartesian to cylindrical coordinates system

    Notes
    -----
    See cylindrical_to_cartesian for more information.

    Parameters
    ----------
    x : float
    y : float
    z : float
    radians : bool
        True if returned azimuth is in radians instead of degrees.

    Returns
    -------
    radius : float
    azimuth : float
    z : float
        Cylindrical coordinate parameters. See cylindrical_to_cartesian for more information.
    """
    cdef double radius, azimuth
    radius = sqrt(x**2 + y**2)
    azimuth = atan2(y, x)
    if not radians:
        azimuth = _rad_to_deg(azimuth)
    return radius, azimuth, z


# ***************
# Python geometry
# ***************

def intersection_start_and_end(subject, reference):
    """ Find start and end indices of the part of the list of points 'subject' that will intersect with 'reference'

    Notes
    -----
    This is a function designed to find assemblies from within bundles of helices.
    The idea here is to take two lists of points, a subject and a reference. We want to find the part of the subject
    (this could be all of it, or none of it) that would intersect the line joining the points in the reference, if
    the line joining the points of subject was the prinicipal axis of a cylinder of infinite radius.
    If subject was a list of ten points joining (0, 0 , 0) and (0, 0, 10), and reference was a list of 1000 points
    joining (10, 0, 7) and (10, 0, 12) this function would return (7, 10).
    If the subject and reference for the above example were swapped around, however, the function would return (0, 599).

    Parameters
    ----------
    subject : list of numpy.array
        A list of points in R^3.
    reference : list of numpy.array
        A list of points in R^3.

    Returns
    -------
    start_index : int
        The index of the first point in subject that is part of the intersection described above.
    end_index : int
        The index of the final point in subject that is part of the intersection described above.

    """

    subject_direction = unit_vector(subject[-1] - subject[0])
    reference_direction = unit_vector(reference[-1] - reference[0])

    # If the directions are obtuse, flip the reference.
    if not is_acute(subject_direction, reference_direction):
        reference = numpy.flipud(reference)

    # Calculate the distances from the start and end of a list to all points in the other list.
    subject_start_to_reference = [distance(subject[0], x) for x in reference]
    subject_end_to_reference = [distance(subject[-1], x) for x in reference]
    reference_start_to_subject = [distance(reference[0], x) for x in subject]
    reference_end_to_subject = [distance(reference[-1], x) for x in subject]

    # If the start of the subject is less remote than the start of the reference then the start_index is 0.
    # Otherwise it is the index of the point in subject that is closest to the start of reference.
    if min(subject_start_to_reference) <= min(reference_start_to_subject):
        start_index = 0
    else:
        start_index = reference_start_to_subject.index(min(reference_start_to_subject))

    # If the end of the subject is less remote than the end of the reference then
    #  the end_index is the index of the final element of subject.
    # Otherwise it is the index of the point in subject that is closest to the end of reference.
    if min(subject_end_to_reference) <= min(reference_end_to_subject):
        end_index = len(subject) - 1
    else:
        end_index = reference_end_to_subject.index(min(reference_end_to_subject))

    return start_index, end_index


def closest_distance(points1, points2):
    """ Minimun distance between two lists of points in R^3.

    Distances are calculated pairwise between points in points1 and points2.
    The minimum of all these distance values is returned.

    Parameters
    ----------
    points1 : list or tuple or numpy.array of triples.
        A number of points in R^3.
    points2: list or tuple or numpy.array of triples.
        A number of points in R^3.

    Returns
    -------
    min_dist : float
        The minimum distance between any point from points1 and any point from point2.
    """
    # Use itertools.product to calculate all pair-wise distances between points1 and points2.
    min_dist = min(distance(x, y) for x, y in itertools.product(points1, points2))
    return min_dist


def points_on_a_circle(n, radius=1, centre=(0, 0), rotation=0):
    """ List of n uniformly distributed (x, y) coordinates on the circumference of a circle.

    Parameters
    ----------
    n : int
        Number of points to return.
    radius : float
    centre : tuple or list or numpy.array
    rotation : float
        Angle in degrees by which all points will be rotated.
        rotation = 0 means that the line from the centre to the first point is parallel to the x-axis.
        rotation > 0 => anti-clockwise rotation

    Returns
    -------
    points : list(2-tuples)
        (x, y) coordinates of uniformly distributed points on the circumference of the circle
    """
    rotation = numpy.deg2rad(rotation)
    thetas = [numpy.divide(i * numpy.pi * 2, n) + rotation for i in range(n)]
    points = [(radius * numpy.cos(theta) + centre[0], radius * numpy.sin(theta) + centre[1]) for theta in thetas]
    return points


def centre_of_mass(points, masses=None):
    """ Find centre of mass of list of points.

    Notes
    -----
    Uniform weighting applied if masses not supplied.

    Parameters
    ----------
    points : list(tuple, list or numpy.array)
    masses : list(float)

    Returns
    -------
    com : numpy.array
        Centre of mass of points.

    Raises
    ------
    ValueError
        if len(masses) != len(points)
    """
    n = len(points)
    if not masses:
        masses = [1 / n] * n
    elif len(masses) != n:
        raise ValueError('Number of points ({0}) must be equal to number of masses supplied ({1})'
                         .format(n, len(masses)))

    points = [numpy.array(p) for p in points]
    total_mass = sum(masses)
    com = sum([m * p for m, p in zip(masses, points)]) / total_mass
    return com


def rmsd(points1, points2):
    """Returns average distance between points in points1 and their counterparts in points2.

    Parameters
    ----------
    points1 : list of [float, float, float]
    points2 : list of [float, float, float]

    Returns
    -------
    rmsd : float

    Raises
    ------
    ValueError : if len(points1) != len(points2)
    """
    if len(points1) != len(points2):
        raise ValueError("points1 ({0}) and points2 ({1}) are not of equal length".format(len(points1), len(points2)))
    rmsd = sum([distance(p1, p2) for p1, p2 in zip(points1, points2)]) / len(points1)
    return rmsd


def find_transformations(s1, e1, s2, e2, radians=False):
    """ Returns geometrical transformations required to align the vector from s1 -> e1 with the vector from s2 -> e2

    Notes
    -----
    Applying the transformations returned to the line from s1 -> e2, will align it with the line from s2 -> e2.
    After transformation, s1 will have been moved onto s2.
    Apply the rotation (given by angle, axis and point) first and then the translation (given by translation).

    Parameters
    ----------
    s1 : [float, float, float]
        3D coordinate. Start of line 1
    e1 : [float, float, float]
        3D coordinate. End of line 1
    s2 : [float, float, float]
        3D coordinate. Start of line 2
    e2 : [float, float, float]
        3D coordinate. End of line 2
    radians : bool
        If True, angle will be returned in radians.

    Returns
    -------
    translation : numpy.array(3)
        3D vector. Translation to be applied to s1 to move it onto s2.
    angle : float
        angle for rotation (default: degrees).
    axis : numpy.array(3)
        axis about which to rotate line 1 to align it with line 2.
    point : numpy.array(3)
        a point lying on the rotation axis.
    """
    s1, e1, s2, e2 = [numpy.array(x) for x in [s1, e1, s2, e2]]
    translation = s2 - s1
    v1 = e1 - s1
    v2 = e2 - s2
    angle = angle_between_vectors(v1, v2, radians=radians)
    if not numpy.isclose(angle, 0):
        axis = numpy.cross(v1, v2)
        # if v1 and v2 parallel or antiparallel
        if numpy.allclose(axis, [0, 0, 0]):
            p = numpy.array([1, 0, 0])
            p_angle = angle_between_vectors(v1, p - v1)
            if (numpy.isclose(p_angle, 0.0)) or (numpy.isclose(p_angle, 180.0)):
                p = numpy.array([0, 1, 0])
                p_angle = angle_between_vectors(v1, p - v1)
                if (numpy.isclose(p_angle, 0.0)) or (numpy.isclose(p_angle, 180.0)):
                    p = numpy.array([0, 0, 1])
            foot = find_foot(s1, e1, p)
            axis = unit_vector(p - foot)
    else:
        # arbitrary axis if angle = 0
        axis = numpy.array([0, 0, 1])
    point = s1
    return translation, angle, axis, point


def find_limits(coordinates):
    cdef float x, y, z
    x, y, z = next(coordinates)
    min_x = max_x = x
    min_y = max_y = y
    min_z = max_z = z
    for x, y, z in coordinates:
        if x < min_x:
            min_x = x
        elif x > max_x:
            max_x = x
        if y < min_y:
            min_y = y
        elif y > max_y:
            max_y = y
        if z < min_z:
            min_z = z
        elif z > max_z:
            max_z = z
    x_limits = (min_x, max_x)
    y_limits = (min_y, max_y)
    z_limits = (min_z, max_z)
    return x_limits, y_limits, z_limits


def gen_sectors(coordinates, box_size):
    """Separates the coordinates into overlapping sectors."""
    box_bounds = []
    coordinates, fl_coords = itertools.tee(coordinates)
    for mn, mx in find_limits(fl_coords):
        md = (mx + mn)/2
        spread = (mx - mn)/2
        n_spread = ((spread//box_size) + 2) * box_size
        box_bounds.append(
            [x+(numpy.ceil(md/box_size)*box_size) for x in numpy.arange(-n_spread, n_spread, box_size)])
    x_bounds = [(x-box_size, x+box_size) for x in box_bounds[0]]
    y_bounds = [(y-box_size, y+box_size) for y in box_bounds[1]]
    z_bounds = [(z-box_size, z+box_size) for z in box_bounds[2]]

    x_sectors = {k:[] for k in x_bounds}
    for coordinate in coordinates:
        for min_max in x_sectors:
            mn, mx = min_max
            if mn < coordinate[0] < mx:
                x_sectors[min_max].append(coordinate)
    for x_sector, coordinates in x_sectors.items():
        y_sectors = {k:[] for k in y_bounds}
        for coordinate in coordinates:
            for min_max in y_bounds:
                mn, mx = min_max
                if mn < coordinate[1] < mx:
                    y_sectors[min_max].append(coordinate)
        x_sectors[x_sector] = y_sectors
    sectors = {}
    for x_sector, y_sectors in x_sectors.items():
        for y_sector, coordinates in y_sectors.items():
            z_sectors = {k:[] for k in z_bounds}
            for coordinate in coordinates:
                for min_max in z_bounds:
                    mn, mx = min_max
                    if mn < coordinate[2] < mx:
                        z_sectors[min_max].append(coordinate)
            for z_sector, coordinates in z_sectors.items():
                if coordinates:
                    sectors[(x_sector, y_sector, z_sector)] = coordinates
    return sectors


# TODO add DualQuaternion class to improve usage of find_transoformation - (screw axis representation of rotation and
# translation. See https://en.wikipedia.org/wiki/Screw_axis and https://en.wikipedia.org/wiki/Dual_quaternion.
class Quaternion:

    def __init__(self, r, i, j, k, is_rotation=False):
        self.r = r
        self.i = i
        self.j = j
        self.k = k
        self.is_rotation = is_rotation

    def __repr__(self):
        return "<Quaternion (r, i, j, k): {0:.5f}, {1:.5f}, {2:.5f}, {3:.5f}>".format(self.r, self.i, self.j, self.k)

    @classmethod
    def real_and_vector(cls, r, vector, is_rotation=False):
        instance = cls(r, *vector, is_rotation)
        return instance

    @classmethod
    def angle_and_axis(cls, angle, axis, radians=False):
        instance = cls(angle, *axis)
        return instance.as_rotation(radians=radians)

    @property
    def vector(self):
        return numpy.array([self.i, self.j, self.k])

    @vector.setter
    def vector(self, value):
        self.i, self.j, self.k = value

    @property
    def as_array(self):
        return numpy.array([self.r, self.i, self.j, self.k])

    @property
    def conjugate(self):
        r = self.r
        vector = -self.vector
        return Quaternion(r, *vector)

    @property
    def norm(self):
        return numpy.sqrt(self.r**2 + self.i**2 + self.j**2 + self.k**2)

    def __add__(self, other):
        if type(other) == Quaternion:
            r = self.r + other.r
            vector = self.vector + other.vector
        elif isinstance(other, Number):
            r = self.r + other
            vector = self.vector
        return Quaternion(r, *vector)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == Quaternion:
            r = self.r - other.r
            vector = self.vector - other.vector
        elif isinstance(other, Number):
            r = self.r - other
            vector = self.vector
        return Quaternion(r, *vector)

    def __invert__(self):
        a = self.r**2 + self.i**2 + self.j**2 + self.k**2
        if a != 0:
            r = self.r / a
            vector = -self.vector / a
            return Quaternion(r, *vector)
        else:
            raise ValueError('Quaternion must be non-zero.')

    def __mul__(self, other):
        if type(other) == Quaternion:
            r = (self.r * other.r) - numpy.dot(self.vector, other.vector)
            vector = (self.r * other.vector) + (other.r * self.vector) + numpy.cross(self.vector, other.vector)
            if self.is_rotation and other.is_rotation:
                is_rotation = True
            else:
                is_rotation = False
        elif isinstance(other, Number):
            r = self.r * other
            vector = self.vector * other
            if other == 1.0 and self.is_rotation:
                is_rotation = True
            else:
                is_rotation = False
        else:
            raise ValueError('Multiplication not defined for type {}'.format(type(other)))
        q = Quaternion(r, *vector)
        q.is_rotation = is_rotation
        return q

    def __rmul__(self, other):
        return self.__mul__(other)

    def __matmul__(self, other):
        return self.apply_conjugation(other)

    def __truediv__(self, other):
        if isinstance(other, Number) and other != 0:
            r = self.r / other
            vector = self.vector / other
        else:
            raise TypeError('Expected type int or float, not {0}'.format(type(other)))
        return Quaternion(r, *vector)

    def distance(self, other):
        return (self - other).norm

    def as_rotation(self, radians=False, force=False):
        if force or (not self.is_rotation):
            if not radians:
                angle = numpy.deg2rad(self.r)
            else:
                angle = self.r
            r = numpy.cos(angle / 2.0)
            vector = numpy.sin(angle / 2.0) * unit_vector(self.vector)
            return Quaternion.real_and_vector(r, vector, is_rotation=True)
        else:
            raise AttributeError('Quaternion is already a rotation, use force=True to override.')

    def apply_conjugation(self, other):
        # if other is a vector
        if (type(other) != Quaternion) and (len(other) == 3):
            if abs(self.norm - 1.0) > 0.0001:
                raise ValueError('If applying conjugation to a vector, quaternion must be unit quaternion')
            conjugated_vector = apply_unit_quaternion(other, self.vector, self.r)
            return conjugated_vector
        # fast way if other has 0 real part
        elif (other.r == 0) and (self.norm == 1.0):
            u = self.vector
            v = other.vector
            s = self.r
            conjugated_vector = apply_unit_quaternion(v, u, s)
            return Quaternion.real_and_vector(r=0, vector=conjugated_vector)
        else:
            return self * other * ~self

    def extract_angle_and_axis(self, radians=False):
        """ Returns the angle and axis of the rotation that the Quaternion represents. """
        half_angle = numpy.arccos(self.r)
        sin_half_angle = numpy.sin(half_angle)
        # zero division problem for trivial angle 0. (But also pi, -pi ??)
        axis = self.vector / sin_half_angle
        angle = 2 * half_angle
        if not radians:
            angle = numpy.rad2deg(angle)
        return angle, axis

    def rotate_vector(self, v, point=None):
        if self.is_rotation:
            if point is None:
                return self @ v
            else:
                return (self @ (v - point)) + point
        else:
            raise AttributeError('Quaternion is not a rotation quaternion.')


##### PARAMETERISED CURVES #######
class Curve(object):
    """ Not intended to be instantiated. This is a base class for other Curve classes to inherit from.

    Notes
    -----
    Any curve must have a point method, which returns the coordinates of a point, given a parameter t.
    """
    def point(self, t):
        raise NotImplementedError


class Axis(Curve):
    """ Axis: a straight line curve in 3-D space, defined by its start and end points (3-D vectors). """
    def __init__(self, start, end):
        self.start = numpy.array(start)
        self.end = numpy.array(end)

    def __repr__(self):
        return "<Axis from {0} to {1}>".format(self.start, self.end)

    def point(self, t):
        return self.start + (t * self.direction_vector)

    @property
    def length(self):
        return numpy.linalg.norm(self.end - self.start)

    @property
    def direction_vector(self):
        return unit_vector(self.end - self.start)

    @property
    def midpoint(self):
        return numpy.divide(numpy.add(self.end, self.start), 2)

    @property
    def unit_tangent(self):
        """ Tangent points along the direction of the line. """
        return self.direction_vector

    @property
    def unit_normal(self):
        """ Normal is perpendicular to tangent. See notes.

        Notes
        -----
        Since the curvature of a straight line is zero, the normal and binormal are not well-defined.
        https://en.wikipedia.org/wiki/Differential_geometry_of_curves#Special_Frenet_vectors_and_generalized_curvatures
        Here we have chosen to use the following convention:
            If the Axis is parallel to the z-axis (i.e. unit_tangent = [0, 0, 1]),
             then the normal is [1, 0, 0] (x-axis) and the binormal is [0, 1, 0] (y-axis).
            If the Axis is not parallel to the z-axis, then we can rotate the z-axis so it aligns with our Axis.
            The normal is the same rotation appled to the x-axis.
            It's helpful here to think of the right-hand rule.
            https://upload.wikimedia.org/wikipedia/commons/thumb/d/d2/Right_hand_rule_cross_product.svg/2000px-Right_hand_rule_cross_product.svg.png
            https://en.wikipedia.org/wiki/Cross_product#Definition
            Using the right-hand rule, the Axis direction vector is your index finger,
            the normal vector is the middle finger and the binormal is the thumb.

        """
        if numpy.allclose(self.unit_tangent, [0, 0, 1]):
            normal_vector = numpy.array([1.0, 0.0, 0.0])
        else:
            translation, angle, axis, point = find_transformations([0, 0, 0], [0, 0, 1], self.start, self.end)
            q = Quaternion.angle_and_axis(angle=angle, axis=axis)
            # unit vector function called to help limit rounding errors.
            normal_vector = unit_vector(q.rotate_vector([1, 0, 0]))
        return normal_vector

    @property
    def unit_binormal(self):
        """ Binormal defined as (Tangent X Normal). """
        tan_vector = numpy.cross(self.unit_tangent, self.unit_normal)
        # Should be unit vector by definition but just in case of rounding errors:
        return unit_vector(tan_vector)

    def distance_to_point(self, point):
        q = find_foot(a=self.start, b=self.end, p=point)
        return distance(q, point)

    def arc_length(self, t, t_ref=0):
        """arc length between point t and point t0

        Notes
        -----
        This is a trivial function in this particular case. Implemented for compatability with Curve class.
        Returns arc length along the curve between the points specified by t and t_ref.
        If t_ref > t then the returned value is negative.

        Parameters
        ----------
        t : float
        t_ref : float
            t and t_ref define points on curve (See self.point()).
        """
        return t - t_ref

    def t_from_arc_length(self, arc_length, t_ref=0):
        """ Returns parameter t that corresponds to a point at distance=arc_length along the curve from the curve
         at t_ref.

        Notes
        -----
        This is a trivial function in this particular case. Implemented for compatability with Curve class.

        Parameters
        ----------
        arc_length : float
            distance along curve
        t_ref : float

        Returns
        -------
        float
            Value of parameter t that corresponds to a point at distance=arc_length along the curve from the point
            defined by parameter t_ref.

        """
        return t_ref + arc_length


    def get_coords(self, n_points=10, spacing=1.52, t_ref=0):
        """ List of evenly spaced points along the curve of specified total length.

        Parameters
        ----------
        n_points : int
            The number of points to be returned.
        spacing : float
            arc-length spacing of points to be returned.
        t_ref : The value of t for specifying the first point in the coordinates list.

        Returns
        -------
        coords : list of numpy.arrays
            3D coordinates of evenly spaced points along the curve.
        """
        length = n_points * spacing
        final_t = self.t_from_arc_length(length, t_ref=t_ref)
        dt = spacing
        # ensure there are at most n_points. # adjust final_t to be sure we'll generate at least n_points
        coords = [self.point(t) for t in numpy.arange(t_ref, final_t + dt, dt)][:n_points]
        return coords


# TODO optimise the methods of this class to use cython math library instead of numpy.
class HelicalCurve(Curve):
    """ HelicalCurve, parameterised by a and b. """
    def __init__(self, a, b, handedness='r', t0=0):
        self.a = a
        self.b = b
        self.t0 = t0
        self.axis_start = numpy.array([0, 0, 0])
        self.axis_end = numpy.array([0, 0, 1])
        if (handedness != 'r') and (handedness != 'l'):
            raise ValueError("Handedness must be either \'r\' (right) or \'l\' (left).")
        self.handedness = handedness

    def __repr__(self):
        if self.handedness == 'r':
            h = 'Right'
        elif self.handedness == 'l':
            h = 'Left'
        return '<HelicalCurve ({0}-handed), a = {1:.3f}, b = {2:.3f} >'.format(h, self.a, self.b)

    @classmethod
    def alpha_and_pitch(cls, alpha, pitch, handedness='r', t0=0, radians=False):
        if not radians:
            alpha = numpy.deg2rad(alpha)
        b = pitch / (2 * numpy.pi)
        a = pitch / (2 * numpy.pi * numpy.tan((numpy.pi / 2) - alpha))
        return cls(a=a, b=b, handedness=handedness, t0=t0)

    @classmethod
    def alpha_and_radius(cls, alpha, radius, handedness='r', t0=0, radians=False):
        if not radians:
            alpha = numpy.deg2rad(alpha)
        a = radius
        b = radius * numpy.tan((numpy.pi / 2) - alpha)
        return cls(a=a, b=b, handedness=handedness, t0=t0)

    @classmethod
    def pitch_and_radius(cls, pitch, radius, handedness='r', t0=0):
        a = radius
        b = pitch / (2 * numpy.pi)
        return cls(a=a, b=b, handedness=handedness, t0=t0)

    @property
    def axis(self):
        return Axis(start=self.axis_start, end=self.axis_end)

    @property
    def radius(self):
        return self.a

    @property
    def pitch(self):
        return 2 * numpy.pi * self.b

    @property
    def alpha(self):
        """ Angle between tangent to helix and the helix axis, in degrees. """
        return numpy.rad2deg((numpy.pi / 2) - numpy.arctan2(self.b, self.a))

    def x(self, t):
        x = self.a * numpy.cos(t)
        if self.handedness == 'l':
            x = -x
        return x

    def y(self, t):
        return self.a * numpy.sin(t)

    def z(self, t):
        return self.b * t

    def point(self, t):
        return numpy.array([self.x(t), self.y(t), self.z(t)])

    def _dx(self, t):
        """ First derivative of x(t) evaluated at t"""
        dx = -self.a * numpy.sin(t)
        if self.handedness == 'l':
            dx = -dx
        return dx

    def _dy(self, t):
        """ First derivative of y(t) evaluated at t"""
        return self.a * numpy.cos(t)

    def _dz(self, t):
        """ First derivative of z(t) evaluated at t"""
        return self.b

    def _ddx(self, t):
        """ Second derivative of x(t) evaluated at t"""
        ddx = -self.a * numpy.cos(t)
        if self.handedness == 'l':
            ddx = -ddx
        return ddx

    def _ddy(self, t):
        """ Second derivative of y(t) evaluated at t"""
        return -self.a * numpy.sin(t)

    def _ddz(self, t):
        """ Second derivative of z(t) evaluated at t"""
        return 0

    def unit_tangent(self, t):
        """ Unit vector pointing along the direction of the curve. """
        tan_vector = numpy.array([self._dx(t), self._dy(t), self._dz(t)])
        return unit_vector(tan_vector)

    def unit_normal(self, t):
        """ Unit vector normal to the curve. Together with unit_tangent, defined the osculating plane at t

        Notes
        -----
        See this article for more information on the tangent, normal and binormal vectors.
        https://en.wikipedia.org/wiki/Differential_geometry_of_curves#Special_Frenet_vectors_and_generalized_curvatures
        """
        normal = numpy.array([self._ddx(t), self._ddy(t), self._ddz(t)])
        return unit_vector(normal)

    def unit_binormal(self, t):
        return numpy.cross(self.unit_tangent(t), self.unit_normal(t))

    def curvature(self):
        return abs(self.a) / (self.a**2 + self.b**2)

    def torsion(self):
        tau = self.b / (self.a**2 + self.b**2)
        if self.handedness == 'l':
            tau = -tau
        return tau

    def arc_length(self, t, t_ref=None):
        """arc length between point t and point t0

        Notes
        -----
        Returns arc length along the curve between the points specified by t and t_ref.
        If t_ref > t then the returned value is negative.

        Parameters
        ----------
        t : float
        t_ref : float
            t and t_ref define points on curve (See self.point()).
        """
        if t_ref is None:
            t_ref = self.t0
        return (t - t_ref) * numpy.sqrt(self.a**2 + self.b**2)

    def t_from_arc_length(self, arc_length, t_ref=None):
        """ Returns parameter t that corresponds to a point at distance=arc_length along the curve from the curve
         at t_ref.

        Parameters
        ----------
        arc_length : float
            distance along curve
        t_ref : float

        Returns
        -------
        float
            Value of parameter t that corresponds to a point at distance=arc_length along the curve from the point
            defined by parameter t_ref.

        """
        if t_ref is None:
            t_ref = self.t0
        return (arc_length / numpy.sqrt(self.a**2 + self.b**2)) + t_ref

    def get_coords(self, n_points=10, spacing=1.52, t_ref=None):
        """ List of evenly spaced points along the curve of specified total length.

        Parameters
        ----------
        n_points : int
            The number of points to be returned.
        spacing : float
            arc-length spacing of points to be returned.
        axis_start : list or tuple or numpy.array of length 3
        axis_end : list or tuple or numpy.array of length 3
            Coordinates returned will be for a helical path winding around an axis specified by axis_end - axis_start.
        t_ref : The value of t for specifying the first point in the coordinates list.

        Returns
        -------
        coords : list of numpy.arrays
            3D coordinates of evenly spaced points along the curve.
        """
        if t_ref is None:
            t_ref = self.t0
        length = n_points * spacing
        final_t = self.t_from_arc_length(length, t_ref=t_ref)
        dt = spacing / numpy.sqrt(self.a**2 + self.b**2)
        # ensure there are at most n_points. # adjust final_t to be sure we'll generate at least n_points
        coords = [self.point(t) for t in numpy.arange(t_ref, final_t + dt, dt)][:n_points]
        # translate and rotate coords if (axis_end - axis_start) does not lie on the z-axis
        if (not numpy.allclose(self.axis.unit_tangent, [0, 0, 1])) and (not numpy.allclose(self.axis.start, [0, 0, 0])):
            translation, rotation_angle, rotation_axis, point = find_transformations([0, 0, 0], [0, 0, 1],
                                                                                     self.axis.start, self.axis.end)
            q = Quaternion.angle_and_axis(angle=rotation_angle, axis=rotation_axis)
            coords = [q.rotate_vector(v=c, point=point) + translation for c in coords]
        return coords


__author__ = 'Jack W. Heal, Christopher W. Wood'
__status__ = 'Development'


