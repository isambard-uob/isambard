# distutils: language = c++
#cython: embedsignature=True

from libcpp.string cimport string
from libcpp.vector cimport vector

from functools import total_ordering
import itertools

from tools.geometry import distance


cdef extern from "dataClasses.hpp" namespace "bude":
    cdef cppclass AtomDataFF:
        AtomDataFF(string, string, double, double, double, double, double, double, double) except +
        string atomName, electrostaticType;
        double radius, hydrophobicPotential, hardness, distNpNp, distNpP, radiusScale, electrostaticPotential


cdef extern from "calculateEnergies.cpp" namespace "bude":
    cdef vector[double] calculatePairEnergy(AtomDataFF *, AtomDataFF *, double)


cdef class PyAtomData:
    cdef AtomDataFF *thisptr
    def __cinit__(self, string atomName, string electrostaticType, double radius, double hydrophobicPotential,
                  double hardness, double distNpNp, double distNpP, double radiusScale, double electrostaticPotential):
        self.thisptr = new AtomDataFF(atomName, electrostaticType, radius, hydrophobicPotential, hardness, distNpNp,
                                      distNpP, radiusScale, electrostaticPotential)
    def __dealloc__(self):
        del self.thisptr


cpdef get_within_ff_cutoff(interaction_pairs, double force_field_cutoff):
    """Finds BUFF interactions for a given ampal object and force field.

    Parameters
    ----------
    interaction_pairs: [Monomer]
        Any AMPAL object with a get_atoms method.
    force_field_cutoff: float
        Minimum interaction distance.

    Returns
    -------
    interaction_set: [(Atom, Atom)]
        All of the atom pairs in range of interacting in BUFF but not within
        covalent bond distance.
    """
    cdef double m_dist
    interactions = []
    for ref_atom_a, ref_atom_b in interaction_pairs:
        m_dist = distance(ref_atom_a._vector, ref_atom_b._vector)
        if m_dist <= force_field_cutoff:
            a_atoms = [atom for atom in ref_atom_a.ampal_parent.atoms.values() if atom._ff_id is not None]
            b_atoms = [atom for atom in ref_atom_b.ampal_parent.atoms.values() if atom._ff_id is not None]
            interactions.extend(itertools.product(a_atoms, b_atoms))
    return interactions


cpdef score_interactions(interactions, ff):
    """Scores a set of interactions.

    Parameters
    ----------
    interactions: [(Atom, Atom)]
        An interaction is defined as a pair of Atoms that have had force field
        parameters assigned to them. This means that they should have a
        'ff_params' tag in their tags.
    ff: BuffForceField
        The force field used for scoring.

    Returns
    -------
    buff_score: BuffScore
        A BuffScore object with information about each of the interactions and
        the atoms involved.
    """
    cdef PyAtomData a_params, b_params
    cdef double dist, ffco
    scores = []
    ffco = ff.distance_cutoff
    ffpsd = ff.parameter_struct_dict
    for a, b in interactions:
        a_params = ffpsd[a._ff_id[0]][a._ff_id[1]]
        b_params = ffpsd[b._ff_id[0]][b._ff_id[1]]
        dist = distance(a._vector, b._vector)
        if dist <= ffco:
            scores.append(calculatePairEnergy(a_params.thisptr, b_params.thisptr, dist))
    buff_score = BuffScore(interactions, scores)
    return buff_score


def find_intra_ampal(ampal, distance_cutoff,
                     no_neighbour_backbone=True, backbone_atoms=('N', 'CA', 'C', 'O')):
    """Finds interactions within an AMPAL object and returns the BUFF score.

    Parameters
    ----------
    ampal: AMPAL object
        Any AMPAL object that inherits from BaseAmpal.
    distance_cutoff: float
        Distance (Å) to find interactions.
    no_neighbour_backbone: bool
        Ignore interactions between neighbouring backbone atoms.
    backbone_atoms: (str)
        A tuple containing labels for the backbone atoms in the polymer/s.

    Returns
    -------
    interactions: [(Atom, Atom)]
        All of the atom pairs in range of interacting in BUFF but not within
        covalent bond distance.

    """
    grp = [mon[mon.reference_atom] for mon in ampal.get_monomers() if hasattr(mon, 'reference_atom')]
    gross_interactions = itertools.combinations(grp, 2)
    interactions = get_within_ff_cutoff(gross_interactions, distance_cutoff)
    if no_neighbour_backbone:
        interactions = list(filter(lambda x: check_if_backbone_neighbours(x, backbone_atoms), interactions))
    return interactions


def check_if_backbone_neighbours(interaction, backbone_atoms):
    atom_1, atom_2 = interaction
    if (atom_1.res_label in backbone_atoms) and (atom_2.res_label in backbone_atoms):
        chain1, res1 = atom_1.unique_id[:2]
        chain2, res2 = atom_2.unique_id[:2]
        if (chain2, int(res2)) == (chain1, int(res1) + 1):
            return False
    return True


def score_intra_ampal(ampal, ff):
    """Returns the BUFF score between all atoms in AMPAL object.

    This is just a convenience function that calls find_intra_ampal
    and score_interactions.

    Parameters
    ----------
    ampal: AMPAL Object
        Any AMPAL object that inherits from BaseAmpal.
    ff: BuffForceField
        The force field used for scoring.

    Returns
    -------
    buff_score: BuffScore
        A BuffScore object with information about each of the interactions and
        the atoms involved.
    """
    interactions = find_intra_ampal(ampal, ff.distance_cutoff)
    buff_score = score_interactions(interactions, ff)
    return buff_score


def find_inter_ampal(ampal_objects, distance_cutoff):
    """Finds interactions between AMPAL objects and returns the BUFF score.

    Parameters
    ----------
    ampal_objects: Assembly or [AMPAL objects]
        Either an assembly or an iterable containing AMPAL objects that inherits from BaseAmpal.
    distance_cutoff: float
        Distance (Å) to find interactions.

    Returns
    -------
    interactions: [(Atom, Atom)]
        All of the atom pairs in range of interacting in BUFF but not within
        covalent bond distance.

    """
    ref_atom_lists = map(
        lambda ampal: [mon[mon.reference_atom] for mon in ampal.get_monomers() if hasattr(mon, 'reference_atom')],
        ampal_objects)
    ampal_pairs = itertools.combinations(ref_atom_lists, 2)
    gross_interactions = itertools.chain(*(itertools.product(*ref_atom_lists) for ref_atom_lists in ampal_pairs))
    interactions = get_within_ff_cutoff(gross_interactions, distance_cutoff)
    return interactions


def score_inter_ampal(ampal_objects, ff):
    """Returns the BUFF score between all the AMPAL objects provided.

    This is just a convenience function that calls find_inter_ampal
    and score_interactions.

    Parameters
    ----------
    ampal_objects: Assembly or [AMPAL objects]
        Either an assembly or an iterable containing AMPAL objects that inherits from BaseAmpal.
    ff: BuffForceField
        The force field used for scoring.

    Returns
    -------
    buff_score: BuffScore
        A BuffScore object with information about each of the interactions and
        the atoms involved.
    """
    interactions = find_inter_ampal(ampal_objects, ff.distance_cutoff)
    buff_score = score_interactions(interactions, ff)
    return buff_score


@total_ordering
class BuffScore(object):
    def __init__(self, interactions, scores):
        """BUFF score containing all the component energies.

        Parameters
        ----------
        interactions: [(Atom, Atom)]
            All of the atom pairs in range of interacting in BUFF but not
            within covalent bond distance.
        scores: [(float, float, float)]
            List of component BUFF scores
        """
        self.inter_scores = list(zip(interactions, scores))
        self.steric = 0
        self.desolvation = 0
        self.charge = 0
        for score in scores:
            self.steric += score[0]
            self.desolvation += score[1]
            self.charge += score[2]
        self.total_energy = sum([sum(x) for x in scores])

    def __getitem__(self, item):
        return self.inter_scores[item]

    def __repr__(self):
        return "<BUFF Score {:.2f}: {:.2f} St | {:.2f} De | {:.2f} Ch>".format(
            self.total_energy, self.steric, self.desolvation, self.charge)

    def __eq__(self, other):
        return self.total_energy == other.total_energy

    def __lt__(self, other):
        return self.total_energy < other.total_energy
