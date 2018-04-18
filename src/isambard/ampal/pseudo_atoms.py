"""Contains AMPAL objects representing pseudo atoms."""

from collections import OrderedDict

from ampal.base_ampal import Atom, Monomer, Polymer, write_pdb
from tools.geometry import distance, radius_of_circumcircle


class PseudoGroup(Polymer):
    """Container for `PseudoMonomer`, inherits from `Polymer`.

    Parameters
    ----------
    monomers : PseudoAtom or [PseudoGroup], optional
        `PseudoMonomer` or list containing `PseudoMonomer` objects to form the
        `PseudoGroup`.
    polymer_id : str, optional
        An ID that the user can use to identify the `PseudoGroup`. This is
        used when generating a pdb file using `PseudoGroup().pdb`.
    ampal_parent : ampal.Assembly, optional
        Reference to `Assembly` containing the `PseudoGroup`.
    sl : int, optional
        The default smoothing level used when calculating the
        backbone primitive.

    Attributes
    ----------
    id : str
        `PseudoGroup` ID
    ampal_parent : ampal.Assembly or None
        Reference to `Assembly` containing the `PseudoGroup`
    molecule_type : str
        A description of the type of `Polymer` i.e. Protein, DNA etc.
    ligands : ampal.LigandGroup
        A `LigandGroup` containing all the `Ligands` associated with this
        `PseudoGroup` chain.
    tags : dict
        A dictionary containing information about this AMPAL object.
        The tags dictionary is used by AMPAL to cache information
        about this object, but is also intended to be used by users
        to store any relevant information they have.
    sl : int
        The default smoothing level used when calculating the
        backbone primitive.

    Raises
    ------
    TypeError
        `Polymer` type objects can only be initialised empty or using
        a `Monomer`.
    """

    def __init__(self, monomers=None, polymer_id=' ', ampal_parent=None, sl=2):
        super().__init__(
            monomers=monomers, polymer_id=polymer_id,
            molecule_type='pseudo_group', ampal_parent=ampal_parent, sl=sl)

    def __repr__(self):
        return '<PseudoGroup chain containing {} {}>'.format(
            len(self._monomers),
            'PseudoMonomer' if len(self._monomers) == 1 else 'PseudoMonomers')


class PseudoMonomer(Monomer):
    """Represents a collection of `PsuedoAtoms`.

    Parameters
    ----------
    pseudo_atoms : OrderedDict, optional
        OrderedDict containing Atoms for the `PsuedoMonomer`. OrderedDict
        is used to maintain the order items were added to the
        dictionary.
    mol_code : str, optional
        One or three letter code that represents the `PsuedoMonomer`.
    monomer_id : str, optional
        String used to identify the `PsuedoMonomer`.
    insertion_code : str, optional
        Insertion code of `PsuedoMonomer`, used if reading from pdb.
    is_hetero : bool, optional
        True if is a hetero atom in pdb. Helps with PDB formatting.
    ampal_parent : ampal.PseudoGroup, optional
        Reference to `PseudoGroup` containing the `PsuedoMonomer`.

    Attributes
    ----------
    mol_code : str
        PDB molecule code that represents the `Nucleotide`.
    insertion_code : str
        Insertion code of `Nucleotide`, used if reading from pdb.
    is_hetero : bool
        True if is a hetero atom in pdb. Helps with PDB formatting.
    states : dict
        Contains an `OrderedDicts` containing atom information for each
        state available for the `Nucleotide`.
    id : str
        String used to identify the `Nucleotide`.
    reference_atom : str
        The key that corresponds to the reference `Atom`. This is used
        by various functions, for example backbone primitives are
        calculated using the `Atom` defined using this key.
    ampal_parent : Polynucleotide or None
        A reference to the `Polynucleotide` containing this `Nucleotide`.
    tags : dict
        A dictionary containing information about this AMPAL object.
        The tags dictionary is used by AMPAL to cache information
        about this object, but is also intended to be used by users
        to store any relevant information they have.

    Raises
    ------
    ValueError
        Raised if `mol_code` is not length 1 or 3.
    """

    def __init__(self, pseudo_atoms=None, mol_code='UNK',
                 monomer_id=' ', insertion_code=' ', ampal_parent=None):
        super(PseudoMonomer, self).__init__(
            atoms=pseudo_atoms, monomer_id=monomer_id,
            ampal_parent=ampal_parent)
        self.mol_code = mol_code
        self.insertion_code = insertion_code
        self.is_hetero = True

    def __repr__(self):
        return '<PseudoMonomer containing {} {}. PseudoMonomer code: {}>'.format(
            len(self.atoms), 'PseudoAtom' if len(self.atoms) == 1 else 'PseudoAtoms', self.mol_code)

    @property
    def pdb(self):
        """Generates a PDB string for the `PseudoMonomer`."""
        pdb_str = write_pdb(
            [self], ' ' if not self.tags['chain_id'] else self.tags['chain_id'])
        return pdb_str


class PseudoAtom(Atom):
    """Object containing 3D coordinates and name.

    Notes
    -----
    Used to represent pseudo atoms (e.g. centre_of_mass) in ISAMBARD.

    Parameters
    ----------
    coordinates : 3D Vector (tuple, list, numpy.array)
        Position of `PseudoAtom` in 3D space.
    element : str
        Element of `PseudoAtom`.
    atom_id : str
        Identifier for `PseudoAtom`, usually a number.
    res_label : str, optional
        Label used in `Monomer` to refer to the `PseudoAtom` type i.e.
        "CA" or "OD1".
    occupancy : float, optional
        The occupancy of the `PseudoAtom`.
    bfactor : float, optional
        The bfactor of the `PseudoAtom`.
    charge : str, optional
        The point charge of the `PseudoAtom`.
    state : str
        The state of this `PseudoAtom`. Used to identify `PseudoAtoms`
        with a number of conformations.
    ampal_parent : ampal.Monomer, optional
       A reference to the `Monomer` containing this `PseudoAtom`.

    Attributes
    ----------
    id : str
        Identifier for `PseudoAtom`, usually a number.
    res_label : str
        Label used in `PseudoGroup` to refer to the `Atom` type i.e. "CA" or "OD1".
    element : str
        Element of `Atom`.
    ampal_parent : ampal.PseudoAtom
       A reference to the `PseudoGroup` containing this `PseudoAtom`.
        number of conformations.
    tags : dict
        A dictionary containing information about this AMPAL object.
        The tags dictionary is used by AMPAL to cache information
        about this object, but is also intended to be used by users
        to store any relevant information they have.
    """

    def __init__(self, coordinates, name='', occupancy=1.0, bfactor=1.0,
                 charge=' ', ampal_parent=None):
        super().__init__(coordinates, element='C', atom_id=' ',
                         occupancy=occupancy, bfactor=bfactor,
                         charge=charge, state='A', ampal_parent=ampal_parent)
        self.name = name

    def __repr__(self):
        return ("<PseudoAtom. Name: {}. Coordinates: "
                "({:.3f}, {:.3f}, {:.3f})>".format(
                    self.name, self.x, self.y, self.z))


class Primitive(PseudoGroup):
    """A backbone path composed of `PseudoAtoms`.

    Parameters
    ----------
    pseudo_atoms : OrderedDict, optional
        OrderedDict containing Atoms for the `PsuedoMonomer`. OrderedDict
        is used to maintain the order items were added to the
        dictionary.
    mol_code : str, optional
        One or three letter code that represents the `PsuedoMonomer`.
    monomer_id : str, optional
        String used to identify the `PsuedoMonomer`.
    insertion_code : str, optional
        Insertion code of `PsuedoMonomer`, used if reading from pdb.
    is_hetero : bool, optional
        True if is a hetero atom in pdb. Helps with PDB formatting.
    ampal_parent : ampal.PseudoGroup, optional
        Reference to `PseudoGroup` containing the `PsuedoMonomer`.

    Attributes
    ----------
    mol_code : str
        PDB molecule code that represents the `Nucleotide`.
    insertion_code : str
        Insertion code of `Nucleotide`, used if reading from pdb.
    is_hetero : bool
        True if is a hetero atom in pdb. Helps with PDB formatting.
    states : dict
        Contains an `OrderedDicts` containing atom information for each
        state available for the `Nucleotide`.
    id : str
        String used to identify the `Nucleotide`.
    reference_atom : str
        The key that corresponds to the reference `Atom`. This is used
        by various functions, for example backbone primitives are
        calculated using the `Atom` defined using this key.
    ampal_parent : Polynucleotide or None
        A reference to the `Polynucleotide` containing this `Nucleotide`.
    tags : dict
        A dictionary containing information about this AMPAL object.
        The tags dictionary is used by AMPAL to cache information
        about this object, but is also intended to be used by users
        to store any relevant information they have.

    Raises
    ------
    ValueError
        Raised if `mol_code` is not length 1 or 3.
    """

    def __init__(self, monomers=None, polymer_id=' ', ampal_parent=None, sl=2):
        super().__init__(
            monomers=monomers, polymer_id=polymer_id,
            ampal_parent=ampal_parent, sl=sl)

    def __repr__(self):
        return '<Primitive chain containing {} {}>'.format(
            len(self._monomers),
            'PseudoMonomer' if len(self._monomers) == 1 else 'PseudoMonomers')

    @classmethod
    def from_coordinates(cls, coordinates):
        """Creates a `Primitive` from a list of coordinates."""
        prim = cls()
        for coord in coordinates:
            pm = PseudoMonomer(ampal_parent=prim)
            pa = PseudoAtom(coord, ampal_parent=pm)
            pm.atoms = OrderedDict([('CA', pa)])
            prim.append(pm)
        prim.relabel_all()
        return prim

    @property
    def coordinates(self):
        """Returns the backbone coordinates for the `Primitive`."""
        return [x._vector for x in self.get_atoms()]

    def rise_per_residue(self):
        """The rise per residue at each point on the Primitive.

        Notes
        -----
        Each element of the returned list is the rise per residue,
        at a point on the Primitive. Element i is the distance
        between primitive[i] and primitive[i + 1]. The final value
        is None.
        """
        rprs = [distance(self[i]['CA'], self[i + 1]['CA'])
                for i in range(len(self) - 1)]
        rprs.append(None)
        return rprs

    def radii_of_curvature(self):
        """The radius of curvature at each point on the Polymer primitive.

        Notes
        -----
        Each element of the returned list is the radius of curvature,
        at a point on the Polymer primitive. Element i is the radius
        of the circumcircle formed from indices [i-1, i, i+1] of the
        primitve. The first and final values are None.
        """
        rocs = []
        for i in range(len(self)):
            if 0 < i < len(self) - 1:
                rocs.append(radius_of_circumcircle(
                    self[i - 1]['CA'], self[i]['CA'], self[i + 1]['CA']))
            else:
                rocs.append(None)
        return rocs


__author__ = 'Jack W. Heal'
