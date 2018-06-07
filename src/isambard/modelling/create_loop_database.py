"""Module for creating loop databases used in ISAMBARD loop modeller."""

import argparse
import glob
import itertools
import multiprocessing
import pathlib
from typing import List, Tuple, Union

from ampal import load_pdb
from ampal.align import align_backbones
from sqlalchemy import create_engine, Column, Float, Integer, String, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from .extract_loop_data import gather_loops_from_pdb


BASE = declarative_base()


class Loop(BASE):
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
            remove_redundant(args.loop_db, args.remove_redundant, args.verbose)
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
            loop_lists = map(fault_tolerant_gather, paths)
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


def fault_tolerant_gather(path: str) -> List[Union[dict, Tuple[str, Exception]]]:
    """Returns the exception rather than raising on failure."""
    try:
        loops = gather_loops_from_pdb(path)
    except Exception as e:
        loops = [(path, e)]
    return loops


def remove_redundant(database_path: str, redundancy_cutoff: float,
                     verbose: bool) -> None:
    loop_db = create_db_session(database_path)
    loop_matches = loop_db.query(Loop)
    loop_sequences: dict = {}
    db_stats = {
        'total_loops': 0,
        'non-redundant': 0,
        'redundant': 0,
    }
    for loop in loop_matches:
        db_stats['total_loops'] += 1
        if loop.sequence in loop_sequences:
            redundant = False
            mobile_ampal = load_pdb(loop.coordinates, path=False)
            for loop_id in loop_sequences[loop.sequence]:
                reference_loop = loop_db.query(Loop).get(loop_id)
                ref_ampal = load_pdb(reference_loop.coordinates,
                                     path=False)
                rmsd = align_backbones(mobile_ampal, ref_ampal,
                                       stop_when=redundancy_cutoff,
                                       verbose=verbose)
                if rmsd <= redundancy_cutoff:
                    redundant = True
                    break
            if not redundant:
                db_stats['non-redundant'] += 1
                loop_sequences[loop.sequence].append(loop.id)
            else:
                if verbose:
                    print(f'\n{loop} is redundant.')
                db_stats['redundant'] += 1
                loop_db.delete(loop)
        else:
            db_stats['non-redundant'] += 1
            loop_sequences[loop.sequence] = [loop.id]
    print(('\n{total_loops} analysed: {non-redundant} non-redundant, {redundant} '
           'redundant.').format(**db_stats))
    loop_db.commit()
    return


def _align_eval(loop, reference):
    return loop.rmsd(reference, backbone=True)


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
              "creating your database to reduce this."
              ),
        metavar='CUT-OFF',
        type=float
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
