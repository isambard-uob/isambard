"""A module for finding and fitting loops regions to a model."""

import copy
import pathlib
from typing import List

# Types are ignored as either no stubs or Cython code
import numpy  # type: ignore
from ampal import geometry, Polypeptide, Residue, load_pdb  # type: ignore
from ampal.align import MMCAlign  # type: ignore

from .create_loop_database import create_db_session, Loop
from .extract_loop_data import calculate_loop_geometry


def find_loops(
        entering_residue: Residue, exiting_residue: Residue, path_to_db: str,
        distance_threshold: float=1.0, angle_threshold: float=10.0,
        dihedral_threshold: float=10.0, loop_type: str='%%', min_length: int=1,
        max_length: int=20, max_resolution: float=3.0) -> List[Loop]:
    """Find loops that fit between the entering and exiting residue.

    Parameters
    ----------
    entering_residue : Residue
        The residue that should enter the loop.
    exiting_residue : Residue
        The residue that should exit the loop.
    path_to_db : str
        Path to a loop database.
    distance_threshold : float, optional
        The variation that is allowed on the distance between the
        entering and exiting residue of the loop.
    angle_threshold : float, optional
        The variation that is allowed on the angle between the
        entering and exiting backbone primitives relative to the
        vector between the entering and exiting residue.
    dihedral_threshold : float, optional
        The variation that is allowed on the dihedral angle between
        the planes formed by average primitive of the 4 entering residues,
        the entering residue primitive and the exiting primitive; and the
        entering primitive, exiting primitive and the average of the
        primitives of the 4 exiting residues.
    loop_type : str, optional
        The DSSP assignment of the secondary structure entering and
        exiting the loop i.e. 'HH' is helix loop helix, 'HE' is helix
        loop beta.
    min_length : int, optional
        The minimum loop length.
    max_length : int, optional
        The maximum loop length.
    max_resolution : float, optional
        The maximum allowed loop resolution i.e. 3.0 will use loops
        with 3.0 A resolution and below.

    Returns
    -------
    loops
        A list of loop records from the loop database.

    Raises
    ------
    ValueError
        Raised if the entering or exiting residues have less than 3
        residues before and after them respectively.

        Raised if the loop database does not exist at the path supplied.
    """
    path = pathlib.Path(path_to_db)
    try:
        assert path.exists() and path.is_file()
    except AssertionError:
        raise ValueError(
            'path_to_db should be a path to a loop database created with the '
            '`mk_loop_db` tool packaged with ISAMBARD. Please check the path '
            'and run again.')
    ent_index = entering_residue.parent._monomers.index(entering_residue)
    exi_index = exiting_residue.parent._monomers.index(exiting_residue)
    try:
        entering_prims = entering_residue.parent.primitive[ent_index-3:ent_index+1]
        exiting_prims = exiting_residue.parent.primitive[exi_index:exi_index+4]
        loop_geometry = calculate_loop_geometry(entering_prims, exiting_prims)
    except IndexError:
        raise ValueError(
            'The entering residue must have at least 3 residues before it and '
            'the exiting residue must have at least 3 residues after it in '
            'order to determine loop geometry.'
        )
    loops = query_loop_database(
        path_to_db, loop_geometry, distance_threshold, angle_threshold,
        dihedral_threshold, loop_type, min_length, max_length, max_resolution)
    return loops


def query_loop_database(
        path_to_db: str, loop_geometry: dict, distance_threshold: float,
        angle_threshold: float, dihedral_threshold: float, loop_type: str,
        min_length: int, max_length: int, max_resolution: float):
    """Queries the loop database for loops with matching geometry."""
    loop_db_session = create_db_session(path_to_db)
    loop_matches = loop_db_session.query(Loop).filter(
        Loop.end_to_end_distance.between(
            loop_geometry['end_to_end_distance'] - distance_threshold,
            loop_geometry['end_to_end_distance'] + distance_threshold
        ),
        Loop.entering_angle.between(
            loop_geometry['entering_angle'] - angle_threshold,
            loop_geometry['entering_angle'] + angle_threshold
        ),
        Loop.exiting_angle.between(
            loop_geometry['exiting_angle'] - angle_threshold,
            loop_geometry['exiting_angle'] + angle_threshold
        ),
        Loop.loop_type.like(f'{loop_type[0]}{loop_type[1]}'),
        Loop.resolution <= max_resolution,
        Loop.length.between(min_length, max_length)
    )
    dihedral_query = Loop.enter_exit_torsion.between(
        loop_geometry['enter_exit_torsion'] - dihedral_threshold,
        loop_geometry['enter_exit_torsion'] + dihedral_threshold
    )
    di_plus_thresh = loop_geometry['enter_exit_torsion'] + \
        dihedral_threshold
    if di_plus_thresh > 180:
        diff = di_plus_thresh - 360
        dihedral_query = dihedral_query | Loop.enter_exit_torsion.between(
            -180, diff)
    di_minus_thresh = loop_geometry['enter_exit_torsion'] - \
        dihedral_threshold
    if di_minus_thresh < -180:
        diff = di_minus_thresh + 360
        dihedral_query = dihedral_query | Loop.enter_exit_torsion.between(
            diff, 180)
    loop_matches = loop_matches.filter(dihedral_query)
    loops = loop_matches.all()
    return loops


def align_loop(entering_residue: Residue, exiting_residue: Residue,
               loop: Loop) -> Polypeptide:
    """Converts a loop database entry to a polypeptide and aligns.""

    Parameters
    ----------
    entering_residue
        The residue that should enter the loop.
    exiting_residue
        The residue that should exit the loop.
    loop : Loop
        Database entry or polypeptide for the loop to be fitted.

    Returns
    -------
    loop_pp : Polypeptide
        An ampal polypeptide for the fitted loop.
    """
    ent_index = entering_residue.parent._monomers.index(entering_residue)
    exi_index = exiting_residue.parent._monomers.index(exiting_residue)
    try:
        entering_pp = entering_residue.parent[ent_index-3:ent_index+1]
        exiting_pp = exiting_residue.parent[exi_index:exi_index+4]
        entering_prims = entering_residue.parent.primitive[ent_index-3:ent_index+1]
        exiting_prims = exiting_residue.parent.primitive[exi_index:exi_index+4]
        loop_geometry = calculate_loop_geometry(entering_prims, exiting_prims)
    except IndexError:
        raise ValueError(
            'The entering residue must have at least 3 residues before it and '
            'the exiting residue must have at least 3 residues after it in '
            'order to determine loop geometry.'
        )
    loop_pp = load_pdb(loop.coordinates, path=False)[0]
    first_prim = loop.get_first_primitive()
    translation, angle, axis, point = geometry.find_transformations(
        loop.get_entering_primitive(),
        loop.get_exiting_primitive(),
        loop_geometry['entering_primitive'],
        loop_geometry['exiting_primitive']
    )
    quaterion = geometry.Quaternion.angle_and_axis(
        angle=angle, axis=axis, radians=False)
    for atom in loop_pp.get_atoms():
        atom._vector = quaterion.rotate_vector(v=atom._vector, point=point)
    first_prim = quaterion.rotate_vector(v=first_prim, point=point)
    loop_pp.translate(translation)
    first_prim += translation
    rot_dihedral = geometry.dihedral(
        loop_geometry['first_primitive'],
        loop_geometry['entering_primitive'],
        loop_geometry['exiting_primitive'],
        first_prim
    )
    loop_pp.rotate(
        rot_dihedral,
        loop_geometry['exiting_primitive'] -
        loop_geometry['entering_primitive'],
        loop_geometry['entering_primitive'],
    )
    entering_rmsd = loop_pp[:4].rmsd(entering_pp, backbone=True)
    exiting_rmsd = loop_pp[-4:].rmsd(exiting_pp, backbone=True)
    loop_pp.tags['loop_data'] = loop
    loop_pp.tags['loop_align'] = numpy.mean([entering_rmsd, exiting_rmsd])
    return loop_pp


def fit_loop(entering_residue: Residue, exiting_residue: Residue,
             loop_pp: Polypeptide, fitting_rounds: int=1000,
             temp: float=298.15, max_dist_step: float=1.0,
             max_angle_step: float=10.0, verbose=True) -> Polypeptide:
    """Fits a polypeptide between two residues.

    Parameters
    ----------
    entering_residue
        The residue that should enter the loop.
    exiting_residue
        The residue that should exit the loop.
    loop : Loop or Polypeptide
        Database entry or polypeptide for the loop to be fitted.
    fitting_rounds : int, optional
        Number of Monte Carlo moves to be performed during fitting.
    temp : float, optional, optional
        Temperature used during fitting process.
    max_angle : float, optional
        The maximum variation in rotation that can moved per
        step.
    max_distance : float, optional
        The maximum distance the can be moved per step.
    verbose : bool, optional
        Prints verbose output if true.

    Returns
    -------
    loop_model : Polypeptide
        An ampal polypeptide for the fitted loop.
    """
    ent_index = entering_residue.parent._monomers.index(entering_residue)
    exi_index = exiting_residue.parent._monomers.index(exiting_residue)
    entering_pp = entering_residue.parent[ent_index-3:ent_index+1]
    exiting_pp = exiting_residue.parent[exi_index:exi_index+4]
    fitter = MMCAlign(_loop_eval, [entering_pp, exiting_pp], loop_pp)
    fitter.start_optimisation(fitting_rounds, max_angle_step, max_dist_step,
                              temp=temp, verbose=verbose)
    loop_model = fitter.best_model
    loop_model.tags['loop_fit_quality'] = fitter.best_energy
    return loop_model


def _loop_eval(working_model, entering_pp, exiting_pp):
    """Function used to evaluate the loop after every move."""
    entering_rmsd = working_model[:4].rmsd(entering_pp, backbone=True)
    exiting_rmsd = working_model[-4:].rmsd(exiting_pp, backbone=True)
    return numpy.mean([entering_rmsd, exiting_rmsd])


def merge_loop(entering_residue: Residue, exiting_residue: Residue,
               loop: Polypeptide, use_loop_flanking: bool=False) -> Polypeptide:
    """Creates a polypeptide by merging a loop with 2 polypeptides.

    Parameters
    ----------
    loop : Polypeptide
        Loop as a polypeptide, usually this will be produced by the
        `fit_loop` function.
    entering_residue : Residue
        The residue that should enter the loop. All of the polypeptide
        until this point will be merged.
    exiting_residue : Residue
        The residue that should exit the loop. All of the polypeptide
        after this point will be merged.
    use_loop_flanking : bool, optional
        If true, the entering and exiting regions from the loop will be
        used in place of the entering and exiting regions from the
        references.

    Returns
    -------
    merged_polypeptide : Polypeptide
        The merged polypeptide model.
    """
    ent_index = entering_residue.parent._monomers.index(entering_residue)
    exi_index = exiting_residue.parent._monomers.index(exiting_residue)
    if not use_loop_flanking:
        entering_pp = copy.deepcopy(entering_residue.parent[:ent_index+1])
        exiting_pp = copy.deepcopy(exiting_residue.parent[exi_index:])
        loop_pp = copy.deepcopy(loop[4:-4])
    else:
        entering_pp = copy.deepcopy(entering_residue.parent[:ent_index-3])
        exiting_pp = copy.deepcopy(exiting_residue.parent[exi_index+4:])
        loop_pp = copy.deepcopy(loop)
    for residue in loop_pp:
        residue.tags['merged_loop'] = True
    merged_polypeptide = entering_pp + loop_pp + exiting_pp
    merged_polypeptide.relabel_all()
    merged_polypeptide.tags = {**loop.tags, **entering_pp.tags,
                               **exiting_pp.tags}
    return merged_polypeptide


__author__ = "Christopher W. Wood"
