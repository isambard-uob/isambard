"""A module for finding and fitting loops regions to a model."""

import copy
import math
import pathlib
import random
import sys
from typing import List, Tuple, Union

from ampal import geometry, Assembly, Polypeptide, Residue, load_pdb
import numpy

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


def fit_loop(entering_residue: Residue, exiting_residue: Residue,
             loop: Union[Loop, Polypeptide], fitting_rounds: int=1000,
             temp: float=298.15, max_dist_step: float=1.0,
             max_angle_step: float=10.0) -> Polypeptide:
    """Takes a loop database entry and fits it between two residues.

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
    temp : float, optional
        Temperature used during fitting process.
    max_angle : float
        The maximum variation in rotation that can moved per
        step.
    max_distance : float
        The maximum distance the can be moved per step.

    Returns
    -------
    loop_model : Polypeptide
        An ampal polypeptide for the fitted loop.
    """
    ent_index = entering_residue.parent._monomers.index(entering_residue)
    exi_index = exiting_residue.parent._monomers.index(exiting_residue)
    entering_pp = entering_residue.parent[ent_index-3:ent_index+1]
    exiting_pp = exiting_residue.parent[exi_index:exi_index+4]
    if isinstance(loop, Loop):
        loop_pp = load_pdb(loop.coordinates, path=False)[0]
    elif isinstance(loop, Polypeptide):
        loop_pp = loop
    else:
        raise TypeError('`loop` must either be a Loop (database entry) or a '
                        'Polypeptide (ampal object).')
    # Initial alignment
    translation, angle, axis, point = geometry.find_transformations(
        loop_pp[:4].centre_of_mass, loop_pp[-4:].centre_of_mass,
        entering_pp.centre_of_mass, exiting_pp.centre_of_mass
    )
    loop_pp.rotate(angle, axis, point)
    loop_pp.translate(translation)
    fitter = MMCFitter(loop_pp, entering_pp, exiting_pp)
    fitter.start_optimisation(fitting_rounds, max_angle_step, max_dist_step,
                              temp=temp)
    loop_model = fitter.best_model
    loop_model.tags['loop_data'] = loop
    loop_model.tags['loop_fit_quality'] = fitter.best_energy
    return loop_model


def merge_loop(loop: Polypeptide, entering_residue: Residue,
               exiting_residue: Residue) -> Polypeptide:
    """Creates a polypeptide by merging a loop with 2 polypeptides.

    Parameters
    ----------
    loop : Polypeptide
        Loop as a polypeptide, usually this will be produced by the
        `fit_loop` function.
    entering_residue
        The residue that should enter the loop. All of the polypeptide
        until this point will be merged.
    exiting_residue
        The residue that should exit the loop. All of the polypeptide
        after this point will be merged.

    Returns
    -------
    merged_polypeptide : Polypeptide
        The merged polypeptide model.
    """
    ent_index = entering_residue.parent._monomers.index(entering_residue)
    exi_index = exiting_residue.parent._monomers.index(exiting_residue)
    entering_pp = copy.deepcopy(entering_residue.parent[:ent_index+1])
    exiting_pp = copy.deepcopy(exiting_residue.parent[exi_index:])
    loop_pp = copy.deepcopy(loop[4:-4])
    for residue in loop_pp:
        residue.tags['merged_loop'] = True
    merged_polypeptide = entering_pp + loop_pp + exiting_pp
    merged_polypeptide.tags = {**loop.tags, **entering_pp.tags,
                               **exiting_pp.tags}
    return merged_polypeptide


class MMCFitter:
    """A loop fitting protocol that uses Metropolis Monte Carlo.

    Parameters
    ----------
    loop_pp : Polypeptide
        An ampal polypeptide containing the model of the loop.
    entering_pp : Polypeptide
        Polypeptide for the entering region of the loop (the 4
        residues immediately before the loop).
    exiting_pp : Polypeptide
        Polypeptide for the exiting region of the loop (the 4
        residues immediately after the loop.
    """

    def __init__(self, loop_pp: Polypeptide, entering_pp: Polypeptide,
                 exiting_pp: Polypeptide) -> None:
        self.current_energy = None
        self.best_energy = None
        self.best_model = None
        self.loop_pp = loop_pp
        self.entering_pp = entering_pp
        self.exiting_pp = exiting_pp

    def _eval_function(self, working_model):
        """Function used to evaluate the loop after every move."""
        entering_rmsd = working_model[:4].rmsd(self.entering_pp, backbone=True)
        exiting_rmsd = working_model[-4:].rmsd(self.exiting_pp, backbone=True)
        return entering_rmsd + exiting_rmsd

    def start_optimisation(self, rounds: int, max_angle: float,
                           max_distance: float, temp: float=298.15):
        """Starts the loop fitting protocol.

        Parameters
        ----------
        rounds : int
            The number of Monte Carlo moves to be evaluated.
        max_angle : float
            The maximum variation in rotation that can moved per
            step.
        max_distance : float
            The maximum distance the can be moved per step.
        temp : float, optional
            Temperature used during fitting process.
        """
        self._generate_initial_score()
        self._mmc_loop(rounds, max_angle, max_distance, temp=temp)
        return

    def _generate_initial_score(self):
        """Runs the evaluation function for the initial pose."""
        self.current_energy = self._eval_function(self.loop_pp)
        self.best_energy = copy.deepcopy(self.current_energy)
        self.best_model = copy.deepcopy(self.loop_pp)
        return

    def _mmc_loop(self, rounds, max_angle, max_distance,
                  temp=298.15, verbose=True):
        """The main Metropolis Monte Carlo loop."""
        current_round = 0
        while current_round < rounds:
            working_model = copy.deepcopy(self.loop_pp)
            random_vector = geometry.unit_vector(
                numpy.random.uniform(-1, 1, size=3))
            mode = random.choice(['rotate', 'translate'])
            if mode == 'rotate':
                random_angle = numpy.random.rand() * max_angle
                working_model.rotate(random_angle, random_vector,
                                     working_model.centre_of_mass)
            else:
                random_translation = random_vector * (numpy.random.rand() *
                                                      max_distance)
                working_model.translate(random_translation)
            proposed_energy = self._eval_function(working_model)
            move_accepted = self.check_move(proposed_energy,
                                            self.current_energy, t=temp)
            if move_accepted:
                self.current_energy = proposed_energy
                if self.current_energy < self.best_energy:
                    self.loop_pp = working_model
                    self.best_energy = copy.deepcopy(self.current_energy)
                    self.best_model = copy.deepcopy(working_model)
            if verbose:
                sys.stdout.write(
                    '\rRound: {}, Current RMSD: {}, Proposed RMSD: {} '
                    '(best {}), {}.       '
                    .format(current_round, self.float_f(self.current_energy),
                            self.float_f(proposed_energy), self.float_f(
                                self.best_energy),
                            "ACCEPTED" if move_accepted else "DECLINED")
                )
                sys.stdout.flush()
            current_round += 1
        return

    @staticmethod
    def float_f(f):
        """Formats a float for printing to std out."""
        return '{:.2f}'.format(f)

    @staticmethod
    def check_move(new, old, t):
        """Determines if a model will be accepted."""
        if (t <= 0) or numpy.isclose(t, 0.0):
            return False
        K_BOLTZ = 1.9872041E-003  # kcal/mol.K
        if new < old:
            return True
        else:
            move_prob = math.exp(-(new - old) / (K_BOLTZ * t))
            if move_prob > random.uniform(0, 1):
                return True
        return False


__author__ = "Christopher W. Wood"
