# Reference AMPAL Objects

This folder contains pickled AMPAL objects.

## Sub Folders

* `non_canonical_amino_acids`
    * Contains reference structures for non-canonical amino acids.
    * AMPAL objects in here should be pickled `ampal.Monomer` objects.
    * If the reference AMPAL object is taken from a natural structure you should follow the recommended naming convention: `{name}_ref_{pdb_file_name}_{chain_index}_{residue_index}.pickle`. For example: `hydroxyproline_ref_1bkv_0_6.pickle`