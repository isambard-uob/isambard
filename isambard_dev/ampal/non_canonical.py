"""Contains tools for working with non-canonical animo acids."""

import pathlib
import pickle

import numpy

from tools.geometry import dihedral, find_transformations


FILE_PATH = pathlib.Path(__file__).parent
REF_PATH = FILE_PATH / 'reference_ampals' / 'non_canonical_amino_acids'


def convert_pro_to_hyp(pro):
    """Converts a pro residue to a hydroxypro residue.

    All metadata associated with the original pro will be lost i.e. tags.
    As a consequence, it is advisable to relabel all atoms in the structure
    in order to make them contiguous.

    Parameters
    ----------
    pro: ampal.Residue
        The proline residue to be mutated to hydroxyproline.
    """
    with open(REF_PATH / 'hydroxyproline_ref_1bkv_0_6.pickle', 'rb') as inf:
        hyp_ref = pickle.load(inf)
    align_nab(hyp_ref, pro)
    to_remove = ['CB', 'CG', 'CD']
    for (label, atom) in pro.atoms.items():
        if atom.element == 'H':
            to_remove.append(label)
    for label in to_remove:
        del pro.atoms[label]
    for key, val in hyp_ref.atoms.items():
        if key not in pro.atoms.keys():
            pro.atoms[key] = val
    pro.mol_code = 'HYP'
    pro.mol_letter = 'X'
    pro.is_hetero = True
    pro.tags = {}
    pro.states = {'A': pro.atoms}
    pro.active_state = 'A'
    for atom in pro.get_atoms():
        atom.ampal_parent = pro
        atom.tags = {'bfactor': 1.0, 'charge': ' ',
                     'occupancy': 1.0, 'state': 'A'}
    return


def align_nab(tar, ref):
    """Aligns the N-CA and CA-CB vector of the target monomer.

    Parameters
    ----------
    tar: ampal.Residue
        The residue that will be aligned to the reference.
    ref: ampal.Residue
        The reference residue for the alignment.
    """
    rot_trans_1 = find_transformations(
        tar['N'].array, tar['CA'].array, ref['N'].array, ref['CA'].array)
    apply_trans_rot(tar, *rot_trans_1)
    rot_ang_ca_cb = dihedral(tar['CB'], ref['CA'], ref['N'], ref['CB'])
    tar.rotate(rot_ang_ca_cb, ref['N'].array - ref['CA'].array, ref['N'].array)
    return


def apply_trans_rot(ampal, translation, angle, axis, point, radians=False):
    """Applies a translation and rotation to an AMPAL object."""
    if not numpy.isclose(angle, 0.0):
        ampal.rotate(angle=angle, axis=axis, point=point, radians=radians)
    ampal.translate(vector=translation)
    return


__author__ = "Christopher W. Wood"
