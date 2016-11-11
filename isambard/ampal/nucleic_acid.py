from ampal.base_ampal import Polymer, Monomer


class Polynucleotide(Polymer):
    def __init__(self, monomers=None, polymer_id=' ', ampal_parent=None, sl=2):
        super().__init__(
            monomers=monomers, polymer_id=polymer_id, molecule_type='nucleic_acid', ampal_parent=ampal_parent, sl=sl)

    @property
    def sequence(self):
        """Returns the sequence of the Polymer as a string.

        Returns
        -------
        sequence : str
            String of the monomer sequence of the polymer.
        """
        seq = [x.mol_code for x in self._monomers]
        return ' '.join(seq)


class Nucleotide(Monomer):
    def __init__(
            self, atoms=None, mol_code='UNK', monomer_id=' ', insertion_code=' ', is_hetero=False, ampal_parent=None):
        super().__init__(atoms, monomer_id, ampal_parent=ampal_parent)
        self.mol_code = mol_code
        self.mol_letter = mol_code[-1]
        self.insertion_code = insertion_code
        self.is_hetero = is_hetero
        self.reference_atom = 'P'
