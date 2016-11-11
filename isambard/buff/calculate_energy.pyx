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


cpdef find_buff_interactions(ampal, ff, internal=False):
    """Finds BUFF interactions for a given ampal object and force field.

    Parameters
    ----------
    ampal: AMPAL object
        Any AMPAL object with a get_atoms method.
    ff: BuffForceField
        The force field is used to derive the minimum interaction distance.

    Returns
    -------
    interaction_set: [(Atom, Atom)]
        All of the atom pairs in range of interacting in BUFF but not within
        covalent bond distance.
    """
    cdef double ffco, m_dist
    ffco = ff.distance_cutoff
    cen_mons = []
    for monomer in ampal.get_monomers():
        if hasattr(monomer, 'reference_atom'):
            cen_mons.append((monomer, monomer.atoms[monomer.reference_atom]))
        # TODO: Add centre of mass as an option here. Needs careful thought so not to miss interactions.
    mon_pairs = []
    interactions = []
    ncaco = ['N', 'CA', 'C', 'O']
    for (monomer_a, ref_atom_a), (monomer_b, ref_atom_b) in itertools.combinations(cen_mons, 2):
        if not internal:
            if monomer_a.ampal_parent != monomer_b.ampal_parent:
                m_dist = distance(ref_atom_a, ref_atom_b)
            else:
                continue
        else:
            m_dist = distance(ref_atom_a, ref_atom_b)
        if m_dist <= ffco:
            a_atoms = [atom for atom in monomer_a.atoms.values() if atom._ff_id is not None]
            b_atoms = [atom for atom in monomer_b.atoms.values() if atom._ff_id is not None]
            if not internal:
                interactions.extend(itertools.product(a_atoms, b_atoms))
            else:
                for interaction in itertools.product(a_atoms, b_atoms):
                    chain1, res1 = interaction[0].unique_id[:2]
                    chain2, res2 = interaction[1].unique_id[:2]
                    if (interaction[0].res_label in ncaco) and (interaction[1].res_label in ncaco):
                        if (chain2, int(res2)) == (chain1, int(res1) + 1):
                            pass
                        else:
                            interactions.append(interaction)
                    else:
                        interactions.append(interaction)
    return interactions

cpdef score_interactions(interactions, ff, threshold = 1.1):
    """Scores a set of interactions.

    Parameters
    ----------
    interactions: [(Atom, Atom)]
        An interaction is defined as a pair of Atoms that have had force field
        parameters assigned to them. This means that they should have a
        'ff_params' tag in their tags.
    ff: BuffForceField
        The force field used for scoring.
    threshold: float
        Cutoff distance for assigning interactions that are covalent bonds.

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


def score_ampal(ampal, ff, threshold=1.1, internal=False):
    """Finds interactions in an AMPAL object and returns the BUFF score.

    Parameters
    ----------
    ampal: AMPAL object
        Any AMPAL object that inherits from BaseAmpal.
    ff: BuffForceField
        Force field for BUFF.
    threshold: float
        Cutoff distance for assigning interactions that are covalent bonds.
    internal: bool
        Measures internal energy if true.

    Returns
    -------
    BUFF_score: BUFFScore
        A BUFFScore object with information about each of the interactions and
        the atoms involved.

    """
    interactions = find_buff_interactions(ampal, ff, internal=internal)
    return score_interactions(interactions, ff, threshold=threshold)


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
