import warnings
import numpy
from collections import Counter

from ampal.pseudo_atoms import Primitive
from tools.geometry import angle_between_vectors, dihedral, distance, find_foot, unit_vector, is_acute
from tools.amino_acids import residue_mwt, water_mass, residue_ext_280, residue_pka, residue_charge,\
    side_chain_dihedrals
from tools.isambard_warnings import NoncanonicalWarning


_nc_warning_str = 'Unnatural amino acid detected, this value may be inaccurate.'


def sequence_molecular_weight(seq):
    """Returns the molecular weight of the polypeptide sequence.

    Units = Daltons"""
    if 'X' in seq:
        warnings.warn(_nc_warning_str, NoncanonicalWarning)
    return sum([residue_mwt[aa] * n for aa, n in Counter(seq).items()]) + water_mass


# TODO How to account for cystine Vs cysteine.
# (when there's a disuplhide bond, ext coefficient changes for the bonding pair to 125.
def sequence_molar_extinction_280(seq):
    """Returns the molar extinction coefficient of the polypeptide sequence at 280 nm.

    Units = M/cm"""
    if 'X' in seq:
        warnings.warn(_nc_warning_str, NoncanonicalWarning)
    return sum([residue_ext_280[aa] * n for aa, n in Counter(seq).items()])


def partial_charge(aa, pH):
    """Calculates the partial charge of the amino acid at a particular pH value."""
    difference = pH - residue_pka[aa]
    if residue_charge[aa] > 0:
        difference *= -1
    ratio = (10 ** difference) / (1 + 10 ** difference)
    return ratio


def sequence_charge(seq, pH=7.4):
    """Calculates the total charge of the input polypeptide sequence."""
    if 'X' in seq:
        warnings.warn(_nc_warning_str, NoncanonicalWarning)
    adj_protein_charge = sum(
        [partial_charge(aa, pH) * residue_charge[aa] * n for aa, n in Counter(seq).items()])
    adj_protein_charge += partial_charge('N-term', pH) * residue_charge['N-term']
    adj_protein_charge += partial_charge('C-term', pH) * residue_charge['C-term']
    return adj_protein_charge


def charge_series(seq, granularity=0.1):
    if 'X' in seq:
        warnings.warn(_nc_warning_str, NoncanonicalWarning)
    ph_range = numpy.arange(1, 13, granularity)
    charge_at_ph = [sequence_charge(seq, ph) for ph in ph_range]
    return ph_range, charge_at_ph


def sequence_isoelectric_point(seq, granularity=0.1):
    """Calculates the isoelectric point of the input polypeptide sequence."""
    if 'X' in seq:
        warnings.warn(_nc_warning_str, NoncanonicalWarning)
    ph_range, charge_at_ph = charge_series(seq, granularity)
    abs_charge_at_ph = [abs(ch) for ch in charge_at_ph]
    pi_index = min(enumerate(abs_charge_at_ph), key=lambda x: x[1])[0]
    return ph_range[pi_index]


def measure_sidechain_torsion_angles(residue, verbose=True):
    """Calculates side-chain dihedral angles for a residue

    Parameters
    ----------
    residue : [Residue]
        ampal.Residue object
    verbose : bool
        If true, tells you when a residue does not have any known dihedral angles to measure

    Returns
    -------
    chi_angles: list of torsion angles, length depends on residue type, in range [-pi, pi]

        [0] = chi1 [if applicable]
        [1] = chi2 [if applicable]
        [2] = chi3 [if applicable]
        [3] = chi4 [if applicable]
    """

    chi_angles = []
    aa = residue.mol_code

    if aa not in side_chain_dihedrals:
        if verbose:
            print("Amino acid {} has no known side-chain dihedral".format(aa))
    else:
        for set_atoms in side_chain_dihedrals[aa]:
            required_for_dihedral = set_atoms[0:4]
            try:
                angle = dihedral(residue[required_for_dihedral[0]]._vector, residue[required_for_dihedral[1]]._vector,
                                 residue[required_for_dihedral[2]]._vector, residue[required_for_dihedral[3]]._vector)
                chi_angles.append(angle)
            except KeyError as k:
                print("{0} atom missing from residue {1} {2} - can't assign dihedral".format(
                    k, residue.mol_code, residue.id))
                chi_angles.append(None)

    return chi_angles


def measure_torsion_angles(residues):
    """Calculates the dihedral angles for a list of backbone atoms.

    Parameters
    ----------
    residues : [Residue]
        List of ampal.Residue objects.

    Returns
    -------
    torsion_angles : list of torsion angle triples.
        One triple for each residue, containing torsion angles in the range [-pi, pi].
            [0] omega
            [1] phi
            [2] psi
        For the first residue, omega and phi are not defined. For the final residue, psi is not defined.

    Raises
    ------
    ValueError
        If the number of input residues is less than 2.
    """
    if len(residues) < 2:
        torsion_angles = [(None, None, None)] * len(residues)
    else:
        torsion_angles = []
        for i in range(len(residues)):
            if i == 0:
                res1 = residues[i]
                res2 = residues[i + 1]
                omega = None
                phi = None
                try:
                    psi = dihedral(res1['N']._vector, res1['CA']._vector, res1['C']._vector, res2['N']._vector)
                except KeyError as k:
                    print("{0} atom missing - can't assign psi".format(k))
                    psi = None
                torsion_angles.append((omega, phi, psi))
            elif i == len(residues) - 1:
                res1 = residues[i - 1]
                res2 = residues[i]
                try:
                    omega = dihedral(res1['CA']._vector, res1['C']._vector, res2['N']._vector, res2['CA']._vector)
                except KeyError as k:
                    print("{0} atom missing - can't assign omega".format(k))
                    omega = None
                try:
                    phi = dihedral(res1['C']._vector, res2['N']._vector, res2['CA']._vector, res2['C']._vector)
                except KeyError as k:
                    print("{0} atom missing - can't assign phi".format(k))
                    phi = None
                psi = None
                torsion_angles.append((omega, phi, psi))
            else:
                res1 = residues[i - 1]
                res2 = residues[i]
                res3 = residues[i + 1]
                try:
                    omega = dihedral(res1['CA']._vector, res1['C']._vector, res2['N']._vector, res2['CA']._vector)
                except KeyError as k:
                    print("{0} atom missing - can't assign omega".format(k))
                    omega = None
                try:
                    phi = dihedral(res1['C']._vector, res2['N']._vector, res2['CA']._vector, res2['C']._vector)
                except KeyError as k:
                    print("{0} atom missing - can't assign phi".format(k))
                    phi = None
                try:
                    psi = dihedral(res2['N']._vector, res2['CA']._vector, res2['C']._vector, res3['N']._vector)
                except KeyError as k:
                    print("{0} atom missing - can't assign psi".format(k))
                    psi = None
                torsion_angles.append((omega, phi, psi))
    return torsion_angles


# TODO: Find a home for this
def cc_to_local_params(pitch, radius, oligo):
    """Returns local parameters for an oligomeric assembly
    
    Parameters
    ----------
    pitch : float
        Pitch of assembly
    radius : float
        Radius of assembly
    oligo : int
        Oligomeric state of assembly
        
    Returns
    -------
    pitchloc : float
        Local pitch of assembly (between 2 adjacent component helices)
    rloc : float
        Local radius of assembly
    alphaloc : float
        Local pitch-angle of assembly
        
    """
    rloc = numpy.sin(numpy.pi / oligo) * radius
    alpha = numpy.arctan((2 * numpy.pi * radius) / pitch)
    alphaloc = numpy.cos((numpy.pi / 2) - ((numpy.pi) / oligo)) * alpha
    pitchloc = (2 * numpy.pi * rloc) / numpy.tan(alphaloc)
    return pitchloc, rloc, numpy.rad2deg(alphaloc)


def residues_per_turn(p):
    """ The number of residues per turn at each Monomer in the Polymer.

    Notes
    -----
    Each element of the returned list is the number of residues per turn, at a point on the Polymer primitive.
    Calculated using the relative positions of the CA atoms and the primitive of the Polymer.
    Element i is the calculated from the dihedral angle using the CA atoms of the Monomers with indices [i, i+1] and
     the corresponding atoms of the primitive.
    The final value is None.

    Parameters
    ----------
    p : ampal.protein.Polypeptide

    Returns
    -------
    rpts : [float]
    """
    cas = p.get_reference_coords()
    prim_cas = p.primitive.coordinates
    dhs = [abs(dihedral(cas[i], prim_cas[i], prim_cas[i + 1], cas[i + 1])) for i in range(len(prim_cas) - 1)]
    rpts = [360.0 / dh for dh in dhs]
    rpts.append(None)
    return rpts


def polymer_to_reference_axis_distances(p, reference_axis, tag=True, reference_axis_name='ref_axis'):
    """ Calculates the distances between the primitive of a Polymer and a reference_axis.

    Notes
    -----
    Distances are calculated between each point of the Polymer primitive and the corresponding point in reference_axis.
    In the special case of the helical barrel, if the Polymer is a helix and the reference_axis represents the centre
    of the barrel, then this function returns the radius of the barrel at each point on the helix primitive.
    The points of the primitive and the reference_axis are run through in the same order, so take care with the relative
    orientation of the reference axis when defining it.

    Parameters
    ----------
    p : A Polymer object.
    reference_axis : list(numpy.array or tuple or list)
        Length of reference_axis must equal length of the Polymer.
        Each element of reference_axis represents a point in R^3.
    tag : bool, optional
        If True, tags the Chain with the reference axis coordinates and each Residue with
        its distance to the ref axis.
        Distances are stored at the Residue level, but refer to distances from the CA atom.
    reference_axis_name : str, optional
        Used to name the keys in tags at Chain and Residue level.

    Returns
    -------
    distances : list(float)

    Raises
    ------
    ValueError
        If the Polymer and the reference_axis have unequal length.
    """
    if not len(p) == len(reference_axis):
        raise ValueError("The reference axis must contain the same number of points as the Polymer primitive.")
    prim_cas = p.primitive.coordinates
    ref_points = reference_axis.coordinates
    distances = [distance(prim_cas[i], ref_points[i]) for i in range(len(prim_cas))]
    if tag:
        p.tags[reference_axis_name] = reference_axis
        monomer_tag_name = 'distance_to_{0}'.format(reference_axis_name)
        for m, d in zip(p._monomers, distances):
            m.tags[monomer_tag_name] = d
    return distances


def crick_angles(p, reference_axis, tag=True, reference_axis_name='ref_axis'):
    """ Calculate the Crick angle for each CA atom in the Polymer, relative to an axis.

    Notes
    -----
    The final value is None, since the angle calculation requires pairs of points on both
     the primitive and reference_axis.

    Parameters
    ----------
    p : A Polymer object.
    reference_axis : list(numpy.array or tuple or list)
        Length of reference_axis must equal length of the Polymer.
        Each element of reference_axis represents a point in R^3.
    tag : bool, optional
        If True, tags the Chain with the reference axis coordinates and each Residue with its Crick angle.
        Crick angles are stored at the Residue level, but are calculated using the CA atom.
    reference_axis_name : str, optional
        Used to name the keys in tags at Chain and Residue level.

    Returns
    -------
    cr_angles : list(float)
        The crick angles in degrees for each CA atom of the Polymer.

    Raises
    ------
    ValueError
        If the Polymer and the reference_axis have unequal length.
    """
    if not len(p) == len(reference_axis):
        raise ValueError("The reference axis must contain the same number of points as the Polymer primitive.")
    prim_cas = p.primitive.coordinates
    p_cas = p.get_reference_coords()
    ref_points = reference_axis.coordinates
    cr_angles = [dihedral(ref_points[i], prim_cas[i], prim_cas[i + 1], p_cas[i])
                 for i in range(len(prim_cas) - 1)]
    cr_angles.append(None)

    if tag:
        p.tags[reference_axis_name] = reference_axis
        monomer_tag_name = 'crick_angle_{0}'.format(reference_axis_name)
        for m, c in zip(p._monomers, cr_angles):
            m.tags[monomer_tag_name] = c

    return cr_angles


def alpha_angles(p, reference_axis, tag=True, reference_axis_name='ref_axis'):
    """ Alpha angle calculated using points on the primitive of helix and axis.

    Notes
    -----
    The final value is None, since the angle calculation requires pairs of points along the primitive and axis.
    This is a generalisation of the calculation used to measure the tilt of a helix in a coiled-coil with respect to the
    central axis of the coiled coil.

    Parameters
    ----------
    p : A Polymer object.
    reference_axis : list(numpy.array or tuple or list)
        Length of reference_axis must equal length of the Polymer.
        Each element of reference_axis represents a point in R^3.
    tag : bool, optional
        If True, tags the Chain with the reference axis coordinates and each Residue with its alpha angle.
        Alpha angles are stored at the Residue level, but are calculated using the CA atom.
    reference_axis_name : str, optional
        Used to name the keys in tags at Chain and Residue level.

    Returns
    -------
    alphas : list of float
        The alpha angle for the Polymer at each point of its primitive, in degrees.

    Raises
    ------
    ValueError
        If the Polymer and the reference_axis have unequal length.

    """
    if not len(p) == len(reference_axis):
        raise ValueError("The reference axis must contain the same number of points as the Polymer primitive.")
    prim_cas = p.primitive.coordinates
    ref_points = reference_axis.coordinates
    alphas = [abs(dihedral(ref_points[i + 1], ref_points[i], prim_cas[i], prim_cas[i + 1]))
              for i in range(len(prim_cas) - 1)]
    alphas.append(None)

    if tag:
        p.tags[reference_axis_name] = reference_axis
        monomer_tag_name = 'alpha_angle_{0}'.format(reference_axis_name)
        for m, a in zip(p._monomers, alphas):
            m.tags[monomer_tag_name] = a

    return alphas


def polypeptide_vector(p, start_index=0, end_index=-1, unit=True):
    """ Vector along the Chain primitive (default is from N-terminus to C-terminus).

    Notes
    -----
    start_index and end_index can be changed to examine smaller sections of the Chain, or reversed to change the
     direction of the vector.

    Parameters
    ----------
    p : A Chain object.
    start_index : int, optional
        Default is 0 (start at the N-terminus of the Chain)
    end_index : int, optional
        Default is -1 (start at the C-terminus of the Chain)
    unit : bool
        If True, the vector returned has a magnitude of 1.

    Returns
    -------
    vector : a numpy.array
        vector has shape (1, 3)
    """
    if len(p) <= 1:
        raise ValueError("Polymer should have length greater than 1. Polymer length = {0}".format(len(p)))
    try:
        prim_cas = p.primitive.coordinates
        direction_vector = prim_cas[end_index] - prim_cas[start_index]
    except ValueError:
        direction_vector = p[end_index]['CA'].array - p[start_index]['CA'].array
    if unit:
        direction_vector = unit_vector(direction_vector)
    return direction_vector


# TODO Change functionality so that a Primitive object is returned.
# (e.g. all CA ALA atoms like with primitives).
def reference_axis_from_chains(chains):
    """Average coordinates from a set of primitives calculated from Chains.

    Parameters
    ----------
    chains : list(Chain)

    Returns
    -------
    reference_axis : numpy.array
        The averaged (x, y, z) coordinates of the primitives for the list of Chains.
        In the case of a coiled coil barrel, this would give the central axis for calculating e.g. Crick angles.

    Raises
    ------
    ValueError :
        If the Chains are not all of the same length.
    """
    if not len(set([len(x) for x in chains])) == 1:
        raise ValueError("All chains must be of the same length")

    # First array in coords is the primitive coordinates of the first chain.
    # The orientation of the first chain orients the reference_axis.
    coords = [numpy.array(chains[0].primitive.coordinates)]
    orient_vector = polypeptide_vector(chains[0])
    # Append the coordinates for the remaining chains, reversing the direction in antiparallel arrangements.
    for i, c in enumerate(chains[1:]):
        if is_acute(polypeptide_vector(c), orient_vector):
            coords.append(numpy.array(c.primitive.coordinates))
        else:
            coords.append(numpy.flipud(numpy.array(c.primitive.coordinates)))

    # Average across the x, y and z coordinates to get the reference_axis coordinates
    reference_axis = numpy.mean(numpy.array(coords), axis=0)
    return Primitive.from_coordinates(reference_axis)


def flip_reference_axis_if_antiparallel(p, reference_axis, start_index=0, end_index=-1):
    """ Reverses reference axis if its direction opposes the direction of the Polymer.

    Notes
    -----
    If the angle between the vector for the Polymer and the vector for the reference_axis is > 90 degrees, then the
    reference axis is reversed. This is useful to run before running polymer_to_reference_axis_distances, crick_angles,
    or alpha_angles.
    For more information on the start and end indices, see chain_vector.

    Parameters
    ----------
    p : A Polymer object.
    reference_axis : list(numpy.array or tuple or list)
        Length of reference_axis must equal length of the Polymer.
        Each element of reference_axis represents a point in R^3.
    start_index : int, optional
        Default is 0 (start at the N-terminus of the Polymer)
    end_index : int, optional
        Default is -1 (start at the C-terminus of the Polymer)

    Returns
    -------
    reference_axis : list(numpy.array or tuple or list)

    """
    p_vector = polypeptide_vector(p, start_index=start_index, end_index=end_index)
    if is_acute(p_vector, reference_axis[end_index] - reference_axis[start_index]):
        reference_axis = numpy.flipud(reference_axis)
    return reference_axis


def make_primitive(cas_coords, window_length=3):
    """ Calculates a running average of cas_coords with a fixed averaging window_length.

    Parameters
    ----------
    cas_coords : list(numpy.array or float or tuple)
        Each element of the list must have length 3.
    window_length : int
        The number of coordinate sets to average each time.

    Returns
    -------
    s_primitive : list(numpy.array)
        Each array has length 3.

    Raises
    ------
    ValueError
        If the length of cas_coords is smaller than the window_length.
    """

    if len(cas_coords) >= window_length:
        primitive = []
        count = 0
        for _ in cas_coords[:-(window_length - 1)]:
            group = cas_coords[count:count + window_length]
            average_x = sum([x[0] for x in group]) / window_length
            average_y = sum([y[1] for y in group]) / window_length
            average_z = sum([z[2] for z in group]) / window_length
            primitive.append(numpy.array([average_x, average_y, average_z]))
            count += 1
    else:
        raise ValueError('A primitive cannot be generated for {0} atoms using a (too large) '
                         'averaging window_length of {1}.'.format(len(cas_coords), window_length))
    return primitive


def make_primitive_smoothed(cas_coords, smoothing_level=2):
    """ Generates smoothed primitive from a list of coordinates.

    Parameters
    ----------
    cas_coords : list(numpy.array or float or tuple)
        Each element of the list must have length 3.
    smoothing_level : int
        Number of times to run the averaging.

    Returns
    -------
    s_primitive : list(numpy.array)
        Each array has length 3.

    Raises
    ------
    ValueError
        If the smoothing level is too great compared to the length of cas_coords.
    """
    try:
        s_primitive = make_primitive(cas_coords)
        for x in range(smoothing_level):
            s_primitive = make_primitive(s_primitive)
    except ValueError:
        raise ValueError('Smoothing level {0} too high, try reducing the number of rounds'
                         ' or give a longer Chain (curent length = {1}).'.format(smoothing_level, len(cas_coords)))
    return s_primitive


def make_primitive_extrapolate_ends(cas_coords, smoothing_level=2):
    """Generates smoothed helix primitives and extrapolates ends lost in the smoothing process.

    Notes
    -----
    From an input list of CA coordinates, the running average is calculated to form a primitive.
    The smoothing_level dictates how many times to calculate the running average.
    A higher smoothing_level generates a 'smoother' primitive - i.e. the points on the primitive more closely fit a
     smooth curve in R^3.
    Each time the smoothing level is increased by 1, a point is lost from either end of the primitive.
    To correct for this, the primitive is extrapolated at the ends to approximate the lost values.
    There is a trade-off then between the smoothness of the primitive and its accuracy at the ends.


    Parameters
    ----------
    cas_coords : list(numpy.array or float or tuple)
        Each element of the list must have length 3.
    smoothing_level : int
        Number of times to run the averaging.

    Returns
    -------
    final_primitive : list(numpy.array)
        Each array has length 3.
    """
    try:
        smoothed_primitive = make_primitive_smoothed(cas_coords, smoothing_level=smoothing_level)
    except ValueError:
        smoothed_primitive = make_primitive_smoothed(cas_coords, smoothing_level=smoothing_level - 1)
    # if returned smoothed primitive is too short, lower the smoothing level and try again.
    if len(smoothed_primitive) < 3:
        smoothed_primitive = make_primitive_smoothed(cas_coords, smoothing_level=smoothing_level - 1)
    final_primitive = []
    for ca in cas_coords:
        prim_dists = [distance(ca, p) for p in smoothed_primitive]
        closest_indices = sorted([x[0] for x in sorted(enumerate(prim_dists), key=lambda k: k[1])[:3]])
        a, b, c = [smoothed_primitive[x] for x in closest_indices]
        ab_foot = find_foot(a, b, ca)
        bc_foot = find_foot(b, c, ca)
        ca_foot = (ab_foot + bc_foot) / 2
        final_primitive.append(ca_foot)
    return final_primitive


__author__ = 'Jack W. Heal, Christopher W. Wood, Gail J. Bartlett, Derek N. Woolfson, Kieran L. Hudson'
