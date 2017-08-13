"""Contains code for analysing chemical interactions in AMPAL objects."""

import itertools
import networkx

from ampal.ampal_databases import element_data
from tools.geometry import distance, gen_sectors

core_components = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'HOH']


class Interaction(object):
    """A container for all types of interaction with donor and acceptor.

    Parameters
    ----------
    a : ampal.Atom
        A member of a pairwise interaction.
    b : ampal.Atom
        A member of a pairwise interaction.
    dist : float
        The distance between `a` and `b`.

    Attributes
    ----------
    a : ampal.Atom
        A member of a pairwise interaction.
    b : ampal.Atom
        A member of a pairwise interaction.
    dist : float
        The distance between `Atom` `a` and `b`.
    """

    def __init__(self, a, b, dist):
        self._a = a
        self._b = b
        self.dist = dist

    def __hash__(self):
        return hash((self._a, self._b))

    def __eq__(self, other):
        return (type(self), self._a, self._b) == (type(other), other._a, other._b)

    def __repr__(self):
        am = self._a.ampal_parent
        ac = am.ampal_parent
        bm = self._b.ampal_parent
        bc = bm.ampal_parent
        return '<Interaction between {} {}{} and {} {}{}>'.format(
            self._a.res_label, am.id, ac.id, self._b.res_label, bm.id, bc.id)


class CovalentBond(Interaction):
    """Defines a covalent bond."""

    @property
    def a(self):
        """One `Atom` involved in the covalent bond."""
        return self._a

    @property
    def b(self):
        """One `Atom` involved in the covalent bond."""
        return self._b

    def __repr__(self):
        am = self._a.ampal_parent
        ac = am.ampal_parent
        bm = self._b.ampal_parent
        bc = bm.ampal_parent
        return '<Covalent bond between {}{} {} {} --- {} {} {}{}>'.format(
            ac.id, am.id, am.mol_code, self._a.res_label, self._b.res_label,
            bm.mol_code, bc.id, bm.id)


class NonCovalentInteraction(Interaction):
    """ A container for all non-covalent interaction.

    Parameters
    ----------
    donor : ampal.Atom
        The donor `Atom` in the interaction.
    acceptor : ampal.Atom
        The acceptor atom in the interaction.
    dist : float
        The distance between `Atom` `a` and `b`.

    Attributes
    ----------
    donor : ampal.Atom
        The donor `Atom` in the interaction.
    acceptor : ampal.Atom
        The acceptor atom in the interaction.
    dist : float
        The distance between `Atom` `a` and `b`.
    """

    def __init__(self, donor, acceptor, dist):
        super().__init__(donor, acceptor, dist)

    @property
    def donor(self):
        """The donor `Atom` in the interaction."""
        return self._a

    @property
    def acceptor(self):
        """The acceptor in the interaction."""
        return self._b

    def __repr__(self):
        return ('<Interaction between {} {}{} (donor) '
                'and {} {}{} (acceptor)>'.format(
                    self.donor.mol_code, self.donor.id, 
                    self.donor.ampal_parent.id, self.acceptor.mol_code,
                    self.acceptor.id, self.acceptor.ampal_parent.id)


class HydrogenBond(NonCovalentInteraction):
    """Defines a hydrogen bond in terms of a donor and an acceptor.

    Parameters
    ----------
    donor : ampal.Atom
        The donor `Atom` in the interaction.
    acceptor : ampal.Atom
        The acceptor atom in the interaction.
    dist : float
        The distance between `Atom` `a` and `b`.
    ang_a : float
        Angle between the acceptor and the interaction vector.
    ang_d : float
        Angle between the donor and the interaction vector.

    Attributes
    ----------
    donor : ampal.Atom
        The donor `Atom` in the interaction.
    acceptor : ampal.Atom
        The acceptor atom in the interaction.
    dist : float
        The distance between `Atom` `a` and `b`.
    ang_a : float
        Angle between the acceptor and the interaction vector.
    ang_d : float
        Angle between the donor and the interaction vector.
    """

    def __init__(self, donor, acceptor, dist, ang_d, ang_a):
        super().__init__(donor, acceptor, dist)
        self.ang_d=ang_d
        self.ang_a=ang_a

    @property
    def donor_monomer(self):
        """The donor `Monomer` in the interaction."""
        return self._a.ampal_parent

    @property
    def acceptor_monomer(self):
        """The acceptor `Monomer` in the interaction."""
        return self._b.ampal_parent

    def __repr__(self):
        dm=self.donor.ampal_parent
        dc=dm.ampal_parent
        am=self.acceptor.ampal_parent
        ac=am.ampal_parent
        return '<Hydrogen Bond between ({}{}) {}-{} ||||| {}-{} ({}{})>'.format(
            dm.id, dc.id, dm.mol_code, self.donor.res_label,
            self.acceptor.res_label, am.mol_code, am.id, ac.id)


def covalent_bonds(atoms, threshold=1.1):
    """Returns all the covalent bonds in a list of `Atom` pairs.

    Notes
    -----
    Uses information `element_data`, which can be accessed directly
    through this module i.e. `isambard.ampal.interactions.element_data`.

    Parameters
    ----------
    atoms : [(`Atom`, `Atom`)]
        List of pairs of `Atoms`.
    threshold : float, optional
        Allows deviation from ideal covalent bond distance to be included.
        For example, a value of 1.1 would allow interactions up to 10% further
        from the ideal distance to be included.
    """
    bonds=[]
    for a, b in atoms:
        bond_distance=(
                element_data[a.element.title()]['atomic radius'] + element_data[
                b.element.title()]['atomic radius']) / 100
        dist=distance(a._vector, b._vector)
        if dist <= bond_distance * threshold:
            bonds.append(CovalentBond(a, b, dist))
    return bonds


def find_covalent_bonds(ampal, max_range=2.2, threshold=1.1, tag=True):
    """Finds all covalent bonds in the AMPAL object.

    Parameters
    ----------
    ampal : AMPAL Object
        Any AMPAL object with a `get_atoms` method.
    max_range : float, optional
        Used to define the sector size, so interactions at longer ranges
        will not be found.
    threshold : float, optional
        Allows deviation from ideal covalent bond distance to be included.
        For example, a value of 1.1 would allow interactions up to 10% further
        from the ideal distance to be included.
    tag : bool, optional
        If `True`, will add the covalent bond to the tags dictionary of
        each `Atom` involved in the interaction under the `covalent_bonds`
        key.
    """
    sectors=gen_sectors(ampal.get_atoms(), max_range * 1.1)
    bonds=[]
    for sector in sectors.values():
        atoms=itertools.combinations(sector, 2)
        bonds.extend(covalent_bonds(atoms, threshold=threshold))
    bond_set=list(set(bonds))
    if tag:
        for bond in bond_set:
            a, b=bond.a, bond.b
            if 'covalent_bonds' not in a.tags:
                a.tags['covalent_bonds']=[b]
            else:
                a.tags['covalent_bonds'].append(b)
            if 'covalent_bonds' not in b.tags:
                b.tags['covalent_bonds']=[a]
            else:
                b.tags['covalent_bonds'].append(a)
    return bond_set


def generate_covalent_bond_graph(covalent_bonds):
    """Generates a graph of the covalent bond network described by the interactions.

    Parameters
    ----------
    covalent_bonds: [CovalentBond]
        List of `CovalentBond`.

    Returns
    -------
    bond_graph: networkx.Graph
        A graph of the covalent bond network.
    """
    bond_graph=networkx.Graph()
    for inter in covalent_bonds:
        bond_graph.add_edge(inter.a, inter.b)
    return bond_graph


def generate_bond_subgraphs_from_break(bond_graph, atom1, atom2):
    """Splits the bond graph between two atoms to producing subgraphs.

    Notes
    -----
    This will not work if there are cycles in the bond graph.

    Parameters
    ----------
    bond_graph: networkx.Graph
        Graph of covalent bond network
    atom1: isambard.ampal.Atom
        First atom in the bond.
    atom2: isambard.ampal.Atom
        Second atom in the bond.

    Returns
    -------
    subgraphs: [networkx.Graph]
        A list of subgraphs generated when a bond is broken in the covalent
        bond network.
    """
    bond_graph.remove_edge(atom1, atom2)
    try:
        subgraphs=list(networkx.connected_component_subgraphs(
            bond_graph, copy=False))
    finally:
        # Add edge
        bond_graph.add_edge(atom1, atom2)
    return subgraphs


__author__='Kieran L. Hudson, Christopher W. Wood, Gail J. Bartlett'
