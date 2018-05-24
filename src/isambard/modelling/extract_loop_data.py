"""Creates data about the properties of loops in a polypeptide."""

from collections import OrderedDict
import pathlib

import numpy
import ampal
from ampal.geometry import angle_between_vectors, dihedral, unit_vector


def make_ss_pattern(regions):
    """Creates a string representing the pattern of secondary structure."""
    ss_types = {
        'H': 'a',
        'E': 'b',
        ' ': 'l',
        'T': 'l',
        'S': 'l',
        'B': 'l',
        'G': '3',
        'I': 'p'
    }
    pattern = [ss_types[region[2]] for region in regions]
    return ''.join(pattern)


def extract_data_for_loops(polypeptide, pdb_code, resolution,
                           create_geometry_path=None):
    """Extracts data from a polypeptide regarding its loops.

    The polypeptide must be tagged with DSSP data to be a valid
    input.

    Parameters
    ----------
    polypeptide : ampal.Polypeptide
        A polypeptide that has been tagged with DSSP data.
    pdb_code : str
        A pdb code for the structure that the PDB code has been
        taken from.
    resolution : float
        The resolution of the structure.
    create_geometry_path : str, optional
        If a path to a directory is provided, a pdb file showing the
        geometry of the loops will be created in that directory.

    Returns
    -------
    loops : [dict]
        Returns a list of dictionaries, one for each loop in the
        structure. Each dictionary contains information about the
        loop. Keys are: pdb_code, resolution, loop_type, start_res,
        end_res, chain, sequence, length, end_to_end_distance,
        entering_angle, exiting_angle, dihedral and coordinates.
    """
    if create_geometry_path:
        create_geometry_path = pathlib.Path(create_geometry_path)
    ss_regions = polypeptide.tags['ss_regions']
    region_groups = [ss_regions[i - 1:i + 2]
                     for i in range(1, len(ss_regions) - 1)]
    ss_pattern = make_ss_pattern(ss_regions)
    pattern_groups = [ss_pattern[i - 1:i + 2]
                      for i in range(1, len(ss_pattern) - 1)]
    # Zip pattern groups with regions filtering for xlx
    loop_regions = [r for (p, r) in zip(pattern_groups, region_groups)
                    if p[1] == 'l']
    pp_primitive = polypeptide.primitive
    assert len(pp_primitive) == len(polypeptide)
    loop_geometry = ampal.Assembly()
    loops = []
    for entering_reg, loop_reg, exiting_reg in loop_regions:
        loop_data, loop_geom = create_loop_dict(
            polypeptide, pp_primitive, entering_reg, loop_reg, exiting_reg,
            create_geometry_path=create_geometry_path)
        loop_data['pdb_code'] = pdb_code
        loop_data['resolution'] = resolution
        loops.append(loop_data)
        if loop_geom:
            loop_geometry.append(loop_geom)
    if create_geometry_path:
        loop_geometry.relabel_all()
        with open(str(create_geometry_path / 'loop_geometry.pdb'), 'w') as outf:
            outf.write(loop_geometry.pdb)
    return loops


def create_loop_dict(polypeptide, pp_primitive,
                     entering_reg, loop_reg, exiting_reg,
                     create_geometry_path=None):
    loop_and_flanking = polypeptide.get_slice_from_res_id(
        str(int(entering_reg[1]) - 3), str(int(exiting_reg[0]) + 3))
    start_index = polypeptide._monomers.index(loop_and_flanking[0])
    end_index = polypeptide._monomers.index(loop_and_flanking[-1])
    loop_region_primitives = pp_primitive[start_index:end_index+1]
    assert len(loop_region_primitives) == len(loop_and_flanking)
    entering_prims = loop_region_primitives[:4]
    exiting_prims = loop_region_primitives[-4:]
    end_to_end_vector = exiting_prims[0]['CA'] - entering_prims[-1]['CA']
    end_to_end_distance = numpy.linalg.norm(end_to_end_vector)
    entering_vector = unit_vector(numpy.mean(
        [entering_prims[i+1]['CA'] - entering_prims[i]['CA']
         for i in range(len(entering_prims) - 1)],
        axis=0))
    exiting_vector = unit_vector(numpy.mean(
        [exiting_prims[i+1]['CA'] - exiting_prims[i]['CA']
         for i in range(len(exiting_prims) - 1)],
        axis=0))
    entering_angle = angle_between_vectors(-entering_vector,
                                           end_to_end_vector)
    exiting_angle = angle_between_vectors(-end_to_end_vector,
                                          exiting_vector)
    enter_exit_torsion = dihedral(
        entering_prims[0]['CA']._vector - entering_vector,
        entering_prims[-1]['CA']._vector,
        exiting_prims[0]['CA']._vector,
        exiting_prims[-1]['CA']._vector + exiting_vector
    )
    loop_data = {
        'loop_type': entering_reg[2] + exiting_reg[2],
        'start_res': loop_reg[0],
        'end_res': loop_reg[1],
        'chain': polypeptide.id,
        'sequence': loop_and_flanking.sequence[4:-4],
        'length': len(loop_and_flanking[4:-4]),
        'end_to_end_distance': end_to_end_distance,
        'entering_angle': entering_angle,
        'exiting_angle': exiting_angle,
        'dihedral': enter_exit_torsion,
        'coordinates': loop_and_flanking.pdb
    }
    if create_geometry_path:
        loop_geom = make_loop_geometry(entering_prims, entering_vector,
                                       exiting_prims, exiting_vector)
    else:
        loop_geom = None
    return loop_data, loop_geom


def make_loop_geometry(entering_prims, entering_vector,
                       exiting_prims, exiting_vector):
    """Creates a Polypeptide that displays the loop geometry."""
    r1 = ampal.Residue(OrderedDict([
        ('CA',
         ampal.Atom(entering_prims[-1]['CA']._vector -
                    (entering_vector*4), element='C'))
    ]))
    r2 = ampal.Residue(OrderedDict([
        ('CA',
         ampal.Atom(entering_prims[-1]['CA']._vector, element='C'))
    ]))
    r3 = ampal.Residue(OrderedDict([
        ('CA',
         ampal.Atom(exiting_prims[0]['CA']._vector, element='C'))
    ]))
    r4 = ampal.Residue(OrderedDict([
        ('CA',
         ampal.Atom(exiting_prims[0]['CA']._vector
                    + (exiting_vector * 4), element='C'))
    ]))
    return ampal.Polypeptide([r1, r2, r3, r4])


__author__ = 'Christopher W. Wood'
