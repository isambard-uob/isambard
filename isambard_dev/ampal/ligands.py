from ampal.base_ampal import Polymer, Monomer


class LigandGroup(Polymer):
    def __init__(self, monomers=None, polymer_id=' ', ampal_parent=None, sl=2):
        super().__init__(
            monomers=monomers, polymer_id=polymer_id, molecule_type='ligands',
            ampal_parent=ampal_parent, sl=sl)

    def __repr__(self):
        return '<Ligands chain containing {} {}>'.format(
            len(self._monomers),
            'Ligand' if len(self._monomers) == 1 else 'Ligands')

    @property
    def categories(self):
        category_dict = {}
        for ligand in self:
            if ligand.category in category_dict:
                category_dict[ligand.category].append(ligand)
            else:
                category_dict[ligand.category] = [ligand]
        return category_dict

    @property
    def category_count(self):
        category_dict = self.categories
        count_dict = {category: len(
            category_dict[category]) for category in category_dict}
        return count_dict


class Ligand(Monomer):
    def __init__(self, mol_code, atoms=None, monomer_id=' ', insertion_code=' ',
                 is_hetero=False, ampal_parent=None):
        """`Monomer` that represents a `Ligand`.

        Notes
        -----
        All `Monomers` that do not have dedicated classes are
        represented using the `Ligand` class.

        Parameters
        ----------
        mol_code : str
            PDB molecule code that represents the monomer.
        atoms : OrderedDict, optional
            OrderedDict containing Atoms for the Monomer. OrderedDict
            is used to maintain the order items were added to the
            dictionary.
        monomer_id : str, optional
            String used to identify the residue.
        insertion_code : str, optional
            Insertion code of monomer, used if reading from pdb.
        is_hetero : bool, optional
            True if is a hetero atom in pdb. Helps with PDB formatting.

        Attributes
        ----------
        mol_code : str
            PDB molecule code that represents the `Ligand`.
        insertion_code : str
            Insertion code of `Ligand`, used if reading from pdb.
        is_hetero : bool
            True if is a hetero atom in pdb. Helps with PDB formatting.
        self.states : dict
            Contains an `OrderedDicts` containing atom information for each
            state available for the `Ligand`.
        id : str
            String used to identify the residue.
        ampal_parent : Polymer or None
            A reference to the `LigandGroup` containing this `Ligand`.
        tags : dict
            A dictionary containing information about this AMPAL object.
            The tags dictionary is used by AMPAL to cache information
            about this object, but is also intended to be used by users
            to store any relevant information they have.

        """
        super(Ligand, self).__init__(
            atoms, monomer_id, ampal_parent=ampal_parent)
        self.mol_code = mol_code
        self.insertion_code = insertion_code
        self.is_hetero = is_hetero

    def __repr__(self):
        return '<Ligand containing {} {}. Ligand code: {}>'.format(
            len(self.atoms), 'Atom' if len(self.atoms) == 1 else 'Atoms',
            self.mol_code)


__author__ = "Christopher W. Wood, Kieran L. Hudson"
