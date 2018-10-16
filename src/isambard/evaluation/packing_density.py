"""Module for evaluating the packing density of a polymer or assembly. """

import ampal
import numpy as np


def tag_packing_density(structure, radius=7):
    """
    Calculates the packing density of each non-hydrogen atom in a polymer
    or assembly.

    An atom's packing density is a measure of the number of atoms within
    its local environment. There are several different methods of
    calculating packing density; we use atomic contact number [1], which
    is the number of non-hydrogen atoms within a specified radius (default
    7 A [1]).

    References
    ----------
    .. [1] Weiss MS (2007) On the interrelationship between atomic
       displacement parameters (ADPs) and coordinates in protein
       structures. *Acta Cryst.* D**63**, 1235-1242.
    """

    if not type(structure).__name__ in ['Polymer', 'Assembly']:
        raise ValueError(
            'Contact order can only be calculated for a polymer or an assembly.'
        )

    atoms_list = [atom for atom in list(structure.get_atoms())
                  if atom.element != 'H']
    atom_coords_array = np.array([atom.array for atom in atoms_list])

    for index, atom in enumerate(atoms_list):
        distances = np.sqrt(np.square(atom_coords_array[:, :] - atom_coords_array[index, :]).sum(axis=1))
        # Subtract 1 to correct for the atom itself being counted
        atom.tags['packing density'] = np.sum(distances < radius) - 1


__author__ = 'Kathryn L. Shelley'
