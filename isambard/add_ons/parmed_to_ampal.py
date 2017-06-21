import parmed
from collections import OrderedDict
import itertools
import tempfile
import os
from copy import deepcopy

from tools.amino_acids import standard_amino_acids
from ampal.ampal_databases import element_data
from ampal.protein import Polypeptide, Residue
from ampal.assembly import Assembly, AmpalContainer
from ampal.base_ampal import Atom
from ampal.ligands import Ligand, LigandGroup
from ampal.nucleic_acid import Nucleotide, Polynucleotide

elements_atomic_numbers = {k: element_data[k]['atomic number'] for k in element_data.keys()}
nucleic_acid_mol_codes = ['A', 'C', 'G', 'U', 'I', 'DA', 'DC', 'DG', 'DT', 'DI']


def convert_cif_to_ampal(cif, path=True, assembly_id=''):
    if not path:
        if type(cif) == str:
            cif = cif.encode()
        try:
            temp_cif = tempfile.NamedTemporaryFile(delete=False)
            temp_cif.write(cif)
            temp_cif.seek(0)
            s = parmed.formats.CIFFile.parse(temp_cif.name, skip_bonds=True)
            temp_cif.close()
        finally:
            os.remove(temp_cif.name)
    else:
        s = parmed.formats.CIFFile.parse(cif, skip_bonds=True)
    return parmed_to_ampal(parmed_structure=s, assembly_id=assembly_id)


def parmed_to_ampal(parmed_structure, assembly_id=''):
    ampal_assembly = Assembly(assembly_id=assembly_id)
    monomer_list = []
    polymer_list = []
    seen_chains = []
    for r in filter(lambda x: not is_ligand(x), parmed_structure.residues):
        if r.chain not in seen_chains:
            seen_chains.append(r.chain)
            if is_protein(r):
                polymer_list.append(Polypeptide(polymer_id=r.chain, ampal_parent=ampal_assembly))
            elif is_nucleic_acid(r):
                polymer_list.append(Polynucleotide(polymer_id=r.chain, ampal_parent=ampal_assembly))
    for r in parmed_structure.residues:
        atoms = OrderedDict([(x.name, parmed_atom_to_ampal_atom(x)) for x in r.atoms])
        polymer = next(filter(lambda x: x.id == r.chain, polymer_list), None)
        if is_protein(r):
            monomer_type = 'p'
        elif is_nucleic_acid(r):
            monomer_type = 'n'
        else:
            monomer_type = 'l'
            if polymer is None:
                polymer = LigandGroup(polymer_id=r.chain, ampal_parent=ampal_assembly)
                polymer_list.append(polymer)
        monomer_dict = dict(mol_code=r.name,
                            monomer_id=r.number,
                            insertion_code=' ' if not r.insertion_code else r.insertion_code,
                            is_hetero=is_hetero(r),
                            ampal_parent=polymer,
                            atoms=atoms)
        if monomer_type == 'p':
            monomer = Residue(**monomer_dict)
        elif monomer_type == 'n':
            monomer = Nucleotide(**monomer_dict)
        else:
            monomer = Ligand(**monomer_dict)
        for x in atoms.values():
            x.ampal_parent = monomer
        monomer_list.append(monomer)
        if isinstance(polymer, LigandGroup):
            polymer.append(monomer)
        else:
            if monomer_type == 'l':
                if polymer.ligands is None:
                    polymer.ligands = LigandGroup()
                polymer.ligands.append(monomer)
                continue
            else:
                polymer.append(monomer)
    ampal_assembly._molecules = polymer_list
    number_of_states = len(parmed_structure.get_coordinates('all'))
    if number_of_states == 1:
        return ampal_assembly
    else:
        ac = AmpalContainer(ampal_objects=[ampal_assembly])
        for i in range(1, number_of_states):
            coords_list = parmed_structure.get_coordinates(i)
            assembly_copy = deepcopy(ampal_assembly)
            for atom in assembly_copy.get_atoms():
                atom_idx = atom.id - 1
                coords = coords_list[atom_idx]
                atom._vector = coords
            ac.append(assembly_copy)
    return ac


def parmed_atom_to_ampal_atom(a, ampal_parent=None):
    c = a.charge
    if c > 0:
        charge = '{0}+'.format(int(c))
    elif c < 0:
        charge = '{0}-'.format(int(c))
    else:
        charge = ''
    ampal_atom = Atom(coordinates=[a.xx, a.xy, a.xz],
                      charge=charge,
                      res_label=a.name,
                      occupancy=a.occupancy,
                      atom_id=a.idx + 1,
                      ampal_parent=ampal_parent,
                      element=get_element(a.element),
                      bfactor=a.bfactor)
    return ampal_atom


def get_element(atomic_number):
    element, at_number = next(filter(lambda x: x[1] == atomic_number, elements_atomic_numbers.items()))
    return element.upper()


def is_protein(parmed_residue):
    names = [x.name for x in parmed_residue.atoms]
    if (all([x in names for x in ['CA', 'N', 'C', 'O']])) or (parmed_residue.name in standard_amino_acids.values()):
        return True
    else:
        return False


def is_nucleic_acid(parmed_residue):
    if parmed_residue.name in nucleic_acid_mol_codes:
        return True
    else:
        return False


def is_hetero(parmed_residue):
    name = parmed_residue.name
    if (name in standard_amino_acids.values()) or (name in nucleic_acid_mol_codes):
        return False
    else:
        return True


def is_ligand(parmed_residue):
    if (not is_protein(parmed_residue)) and (not is_nucleic_acid(parmed_residue)):
        return True
    else:
        return False
