"""Module for creating loop databases used in ISAMBARD loop modeller."""

import argparse
import glob
import itertools
import multiprocessing
import pathlib
from typing import List

# Types are ignored as either no stubs or Cython code
from sqlalchemy import (  # type: ignore
    create_engine, Column, Float, Integer, String, Text)
from sqlalchemy.ext.declarative import declarative_base  # type: ignore
from sqlalchemy.orm import sessionmaker  # type: ignore
from ampal import load_pdb  # type: ignore
from ampal.align import align_backbones  # type: ignore

from .extract_loop_data import gather_loops_from_pdb


BASE = declarative_base()


class Loop(BASE):  # type: ignore
    """Loop database entry."""
    __tablename__ = 'loops'

    id = Column(Integer, primary_key=True)
    pdb_code = Column(String)
    resolution = Column(Float)
    loop_type = Column(String)
    start_res = Column(Integer)
    end_res = Column(Integer)
    chain = Column(String)
    sequence = Column(String)
    length = Column(Integer)
    first_primitive_x = Column(Float)
    first_primitive_y = Column(Float)
    first_primitive_z = Column(Float)
    entering_primitive_x = Column(Float)
    entering_primitive_y = Column(Float)
    entering_primitive_z = Column(Float)
    exiting_primitive_x = Column(Float)
    exiting_primitive_y = Column(Float)
    exiting_primitive_z = Column(Float)
    last_primitive_x = Column(Float)
    last_primitive_y = Column(Float)
    last_primitive_z = Column(Float)
    end_to_end_distance = Column(Float)
    entering_angle = Column(Float)
    exiting_angle = Column(Float)
    enter_exit_torsion = Column(Float)
    coordinates = Column(Text)

    def __repr__(self):
        return "<Loop(pdb='{0}', loop_type='{1}', sequence='{2}')>".format(
            self.pdb_code,
            self.loop_type,
            self.sequence)

    def get_first_primitive(self):
        """Reconstructs the first primitive from its components."""
        primitive = (self.first_primitive_x,
                     self.first_primitive_y,
                     self.first_primitive_z)
        return primitive

    def get_entering_primitive(self):
        """Reconstructs the entering primitive from its components."""
        primitive = (self.entering_primitive_x,
                     self.entering_primitive_y,
                     self.entering_primitive_z)
        return primitive

    def get_exiting_primitive(self):
        """Reconstructs the exiting primitive from its components."""
        primitive = (self.exiting_primitive_x,
                     self.exiting_primitive_y,
                     self.exiting_primitive_z)
        return primitive

    def get_last_primitive(self):
        """Reconstructs the last primitive from its components."""
        primitive = (self.last_primitive_x,
                     self.last_primitive_y,
                     self.last_primitive_z)
        return primitive


def create_db_session(path):
    """Creates a loop database session."""
    engine = create_engine(f'sqlite:///{path}', echo=False)
    BASE.metadata.create_all(engine)
    session = sessionmaker(bind=engine)()
    return session


def main():
    """The main entry point for creating the loop database."""
    args = get_args()
    if args.subparser == 'create':
        data_path = pathlib.Path(args.data_dir)
        data_file_paths = glob.glob(
            str(data_path / '**' / f'*.{args.extension}' if args.recursive
                else data_path / f'*.{args.extension}'),
            recursive=args.recursive)
        process_pdb_files(data_file_paths, args.output, args.processes,
                          args.verbose)

    elif args.subparser == 'fix':
        if args.remove_redundant:
            remove_redundant(args.loop_db, args.remove_redundant,
                             args.processes, args.verbose)
    return


def process_pdb_files(data_file_paths: List[str], output_path: str,
                      processes: int=1, verbose: bool=True) -> None:
    """Extracts the loops from a list of PDB files.

    Parameters
    ----------
    data_file_paths
        A list of paths to files to be processed.
    output_path
        A path (and name) to output the database.
    processes
        The number of processes to be used.
    verbose
        Prints verbose output if true.
    """
    session = create_db_session(output_path)

    db_stats = {
        'total_pdbs': 0,
        'succeeded': 0,
        'failed': 0,
        'total_loops': 0,
    }
    data_chunks = [data_file_paths[i:i + 100]
                   for i in range(0, len(data_file_paths), 100)]
    for paths in data_chunks:
        if processes > 1:
            with multiprocessing.Pool(processes) as pool:
                loop_lists = pool.map(fault_tolerant_gather, paths)
        else:
            loop_lists = list(map(fault_tolerant_gather, paths))
        db_stats['total_pdbs'] += len(paths)
        for loop in itertools.chain(*loop_lists):
            db_stats['total_loops'] += 1
            if isinstance(loop, tuple):
                (path, ex) = loop
                db_stats['failed'] += 1
                if verbose:
                    print(f'Failed to process {path}:\n {str(ex)}\n')
                else:
                    print(f'Failed to process {path}')
                continue
            loop_entry = Loop(**loop)
            session.add(loop_entry)
            db_stats['succeeded'] += 1
        session.commit()
    print('Finished processing {total_pdbs} PDB files:\n'
          'Attempted to add {total_loops} loops to database: {succeeded} '
          'succeeded, {failed} failed.'.format(**db_stats))
    return


def fault_tolerant_gather(path: str) -> list:
    """Returns the exception rather than raising on failure."""
    try:
        loops: list = gather_loops_from_pdb(path)
    except Exception as exception:  # pylint: disable=broad-except
        loops = [(path, exception)]
    return loops


def remove_redundant(database_path: str, redundancy_cutoff: float,
                     processes: int=1, verbose: bool=False) -> None:
    """Removes redundant loops from the loop database."""
    loop_db = create_db_session(database_path)
    loop_matches = loop_db.query(Loop)
    loop_sequences: dict = {}
    for loop in loop_matches:
        if loop.sequence in loop_sequences:
            loop_sequences[loop.sequence].append(loop)
        else:
            loop_sequences[loop.sequence] = [loop]
    multiloops = [loop_sequence for loop_sequence in loop_sequences.values()
                  if len(loop_sequence) > 1]
    if processes > 1:
        with multiprocessing.Pool(processes) as pool:
            redudant_id_lists = pool.map(
                list_redundant,
                [(loop_list, redundancy_cutoff, verbose)
                 for loop_list in multiloops])
    else:
        redudant_id_lists = list(map(
            list_redundant,
            [(loop_list, redundancy_cutoff, verbose)
             for loop_list in multiloops]))
    redundant_ids = list(itertools.chain(*redudant_id_lists))
    redudant_loops = loop_db.query(Loop).filter(Loop.id.in_(redundant_ids))
    redudant_loops.delete(synchronize_session='fetch')
    loop_db.commit()
    print(
        f'\nFound {len(redundant_ids)} redundant loops. Deleted from database.')
    return


def list_redundant(arguments):
    """Creates a list of redundant loops based on 3D alignment."""
    loop_list, redundancy_cutoff, verbose = arguments
    assert len(loop_list) > 1
    non_redundant_loops = [loop_list[0]]
    redundant_loops = []
    for mobile_loop in loop_list[1:]:
        redundant = False
        mobile_ampal = load_pdb(mobile_loop.coordinates, path=False)
        for reference_loop in non_redundant_loops:
            ref_ampal = load_pdb(reference_loop.coordinates, path=False)
            rmsd = align_backbones(mobile_ampal, ref_ampal,
                                   stop_when=redundancy_cutoff,
                                   verbose=False)
            if rmsd <= redundancy_cutoff:
                redundant = True
                break
        if not redundant:
            non_redundant_loops.append(mobile_loop)
        else:
            if verbose:
                print(f'\n{mobile_loop} is redundant.')
            redundant_loops.append(mobile_loop)
    return [x.id for x in redundant_loops]


def get_args():
    """Loads arguments for the main."""
    description = (
        "Commandline tool for creating the loop database for ISAMBARD "
        "loop modelling functionality."
    )
    parser = argparse.ArgumentParser(description=description)
    subparsers = parser.add_subparsers(dest='subparser')
    create = subparsers.add_parser(
        'create', help="Creates a loop database.")
    create.add_argument(
        'data_dir',
        metavar='DATA_DIRECTORY',
        help="Folder containing PDB files used to create the loop database."
    )
    create.add_argument(
        '-o', '--output',
        help=("Path and file name for the database that's created. "
              "If no path is supplied, a file called loops.db will be "
              "created in the current working directory."),
        metavar='OUTPUT_DIRECTORY',
        default='loops.db',
        type=str
    )
    create.add_argument(
        '-e', '--extension',
        help="Extension of structure files to be loaded.",
        metavar='EXTENSION',
        default='pdb',
        type=str
    )
    create.add_argument(
        '-p', '--processes',
        help="Number of processes used to build the loop database.",
        type=int,
        default=1
    )
    create.add_argument(
        '-r', '--recursive',
        help="Will recursively search the DATA_DIRECTORY for PDB files.",
        action="store_true"
    )
    create.add_argument(
        '-v', '--verbose',
        help="Verbose output.",
        action="store_true"
    )
    fix = subparsers.add_parser(
        'fix', help="Fixes a loop database.")
    fix.add_argument(
        'loop_db', metavar="LOOP_DATABASE", help="Path to loop database.")
    fix.add_argument(
        '-r', '--remove-redundant',
        help=("Removes redundant loops from the database with an RMSD above "
              "cut off. Depending on the size of your loop database, this may "
              "take a very long time to finish! It's highly recommended that "
              "you use a non-redundant set of input structure when initially "
              "creating your database to reduce this."),
        metavar='CUT-OFF',
        type=float
    )
    fix.add_argument(
        '-p', '--processes',
        help="Number of processes used to filter redundant loops.",
        type=int,
        default=1
    )
    fix.add_argument(
        '-v', '--verbose',
        help="Verbose output.",
        action="store_true"
    )
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    main()
