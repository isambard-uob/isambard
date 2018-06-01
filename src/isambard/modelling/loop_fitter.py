"""A module for finding and fitting loops regions to a model."""

import pathlib
from typing import List

from ampal import Residue

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
    entering_residue
        The residue that should enter the loop.
    exiting_residue
        The residue that should exit the loop.
    path_to_db
        Path to a loop database.
    distance_threshold, optional
        The variation that is allowed on the distance between the
        entering and exiting residue of the loop.
    angle_threshold, optional
        The variation that is allowed on the angle between the
        entering and exiting backbone primitives relative to the
        vector between the entering and exiting residue.
    dihedral_threshold, optional
        The variation that is allowed on the dihedral angle between
        the planes formed by average primitive of the 4 entering residues,
        the entering residue primitive and the exiting primitive; and the
        entering primitive, exiting primitive and the average of the
        primitives of the 4 exiting residues.
    loop_type, optional
        The DSSP assignment of the secondary structure entering and
        exiting the loop i.e. 'HH' is helix loop helix, 'HE' is helix
        loop beta.
    min_length, optional
        The minimum loop length.
    max_length, optional
        The maximum loop length.
    max_resolution, optional
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
        diff = (di_plus_thresh % 180) - 180
        dihedral_query = dihedral_query | Loop.enter_exit_torsion.between(
            -180, diff)
    di_minus_thresh = loop_geometry['enter_exit_torsion'] - \
        dihedral_threshold
    if di_minus_thresh < -180:
        diff = 180 + (di_minus_thresh % -180)
        dihedral_query = dihedral_query | Loop.enter_exit_torsion.between(
            diff, 180)
    loop_matches.filter(dihedral_query)
    loops = loop_matches.all()
    return loops
