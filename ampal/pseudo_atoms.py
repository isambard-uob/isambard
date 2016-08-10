from collections import OrderedDict

from ampal.base_ampal import Atom, Monomer, Polymer, write_pdb
from tools.geometry import distance, radius_of_circumcircle


class PseudoGroup(Polymer):
    def __init__(self, monomers=None, polymer_id=' ', ampal_parent=None, sl=2):
        super().__init__(
            monomers=monomers, polymer_id=polymer_id, molecule_type='pseudo_group', ampal_parent=ampal_parent, sl=sl)

    def __repr__(self):
        return '<PseudoGroup chain containing {} {}>'.format(
            len(self._monomers), 'PseudoMonomer' if len(self._monomers) == 1 else 'PseudoMonomers')


class PseudoMonomer(Monomer):
    def __init__(self, pseudo_atoms=None, mol_code='UNK',
                 monomer_id=' ', insertion_code=' ', ampal_parent=None):
        """Object containing PseudoAtoms, this is how Monomers are represented in ISAMBARD.

        Parameters
        ----------
        pseudo_atoms : OrderedDict
            OrderedDict containing Atoms for the Monomer. OrderedDict is used to maintain the order items
            were added to the dictionary.
        mol_code : str
            PDB molecule code that represents the monomer.
        monomer_id : str
            String used to identify the residue.
        insertion_code : str
            Insertion code of monomer, used if reading from pdb.
        is_hetero : bool
            True if is a hetero atom in pdb. Helps with PDB formatting.
        """
        super(PseudoMonomer, self).__init__(atoms=pseudo_atoms, monomer_id=monomer_id, ampal_parent=ampal_parent)
        self.mol_code = mol_code
        self.insertion_code = insertion_code
        self.is_hetero = True

    def __repr__(self):
        return '<PseudoMonomer containing {} {}. PseudoMonomer code: {}>'.format(
            len(self.atoms), 'PseudoAtom' if len(self.atoms) == 1 else 'PseudoAtoms', self.mol_code)

    @property
    def pdb(self):
        pdb_str = write_pdb([self], ' ' if not self.tags['chain_id'] else self.tags['chain_id'])
        return pdb_str


class PseudoAtom(Atom):
    """Object containing 3D coordinates and name. Used to represent pseudo atoms (e.g. centre_of_mass) in ISAMBARD."""
    def __init__(self, coordinates, name='', occupancy=1.0, bfactor=1.0, charge=' ', ampal_parent=None):
        super().__init__(coordinates, element='C', atom_id=' ', occupancy=occupancy,
                         bfactor=bfactor, charge=charge, state='A', ampal_parent=ampal_parent)
        self.name = name

    def __repr__(self):
        return "<PseudoAtom. Name: {}. Coordinates: ({:.3f}, {:.3f}, {:.3f})>".format(self.name, self.x, self.y, self.z)


class Primitive(PseudoGroup):
    def __init__(self, monomers=None, polymer_id=' ', ampal_parent=None, sl=2):
        super().__init__(
            monomers=monomers, polymer_id=polymer_id, ampal_parent=ampal_parent, sl=sl)

    def __repr__(self):
        return '<Primitive chain containing {} {}>'.format(
            len(self._monomers), 'PseudoMonomer' if len(self._monomers) == 1 else 'PseudoMonomers')

    @classmethod
    def from_coordinates(cls, coordinates):
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
        return [x._vector for x in self.get_atoms()]

    def rise_per_residue(self):
        """ The rise per residue at each point on the Primitive.

        Notes
        -----
        Each element of the returned list is the rise per residue, at a point on the Primitive.
        Element i is the distance between primitive[i] and primitive[i + 1].
        The final value is None.
        """
        rprs = [distance(self[i]['CA'], self[i + 1]['CA']) for i in range(len(self) - 1)]
        rprs.append(None)
        return rprs

    def radii_of_curvature(self):
        """ The radius of curvature at each point on the Polymer primitive.

        Notes
        -----
        Each element of the returned list is the radius of curvature, at a point on the Polymer primitive.
        Element i is the radius of the circumcircle formed from indices [i-1, i, i+1] of the primitve.
        The first and final values are None.
        """
        rocs = []
        for i in range(len(self)):
            if 0 < i < len(self) - 1:
                rocs.append(radius_of_circumcircle(self[i - 1]['CA'], self[i]['CA'], self[i + 1]['CA']))
            else:
                rocs.append(None)
        return rocs


__author__ = 'Jack W. Heal'


