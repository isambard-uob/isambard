"""Contains AMPAL objects representing nucleic acids."""

from ampal.base_ampal import Polymer, Monomer


class Polynucleotide(Polymer):
    """`Polymer` type object that represents a `Polynucleotide`.

    Parameters
    ----------
    monomers : Nucleotide or [Nucleotide], optional
        `Nucleotide` or list containing `Nucleotide` objects to form the
        `Polynucleotide`.
    polymer_id : str, optional
        An ID that the user can use to identify the `Polynucleotide`. This is
        used when generating a pdb file using `Polynucleotide().pdb`.
    ampal_parent : ampal.Assembly, optional
        Reference to `Assembly` containing the `Polynucleotide`.
    sl : int, optional
        The default smoothing level used when calculating the
        backbone primitive.

    Attributes
    ----------
    id : str
        `Polynucleotide` ID
    ampal_parent : ampal.Assembly or None
        Reference to `Assembly` containing the `Polynucleotide`
    molecule_type : str
        A description of the type of `Polymer` i.e. Protein, DNA etc.
    ligands : ampal.LigandGroup
        A `LigandGroup` containing all the `Ligands` associated with this
        `Polynucleotide` chain.
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
            molecule_type='nucleic_acid', ampal_parent=ampal_parent, sl=sl)

    @property
    def sequence(self):
        """Returns the sequence of the `Polynucleotide` as a string.

        Returns
        -------
        sequence : str
            String of the monomer sequence of the `Polynucleotide`.
        """
        seq = [x.mol_code for x in self._monomers]
        return ' '.join(seq)


class Nucleotide(Monomer):
    """Represents a nucleotide base.

    Parameters
    ----------
    atoms : OrderedDict, optional
        OrderedDict containing Atoms for the `Nucleotide`. OrderedDict
        is used to maintain the order items were added to the
        dictionary.
    mol_code : str, optional
        One or three letter code that represents the `Nucleotide`.
    monomer_id : str, optional
        String used to identify the `Nucleotide`.
    insertion_code : str, optional
        Insertion code of `Nucleotide`, used if reading from pdb.
    is_hetero : bool, optional
        True if is a hetero atom in pdb. Helps with PDB formatting.
    ampal_parent : ampal.Polynucleotide, optional
        Reference to `Polynucleotide` containing the `Nucleotide`.

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

    def __init__(
            self, atoms=None, mol_code='UNK', monomer_id=' ',
            insertion_code=' ', is_hetero=False, ampal_parent=None):
        super().__init__(atoms, monomer_id, ampal_parent=ampal_parent)
        self.mol_code = mol_code
        self.mol_letter = mol_code[-1]
        self.insertion_code = insertion_code
        self.is_hetero = is_hetero
        self.reference_atom = 'P'


__author__ = "Christopher W. Wood"
