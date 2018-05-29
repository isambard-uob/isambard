"""Module for creating loop databases used in ISAMBARD loop modeller."""

import argparse
import glob
import itertools
import multiprocessing
import pathlib
from typing import List, Union

from sqlalchemy import create_engine, Column, Float, Integer, String, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from .extract_loop_data import gather_loops_from_pdb

BASE = declarative_base()


class Loops(BASE):
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
    dihedral = Column(Float)
    coordinates = Column(Text)

    def __repr__(self):
        return "<Loops(pdb='{0}', loop_type='{1}', sequence='{2}) >".format(
            self.pdb,
            self.loop_type,
            self.sequence)


def main():
    """The main entry point for creating the loop database."""
    args = get_args()
    data_path = pathlib.Path(args.data_dir)
    data_file_paths = glob.glob(
        str(data_path / '**' / '*.pdb' if args.recursive
            else data_path / '*.pdb'),
        recursive=args.recursive)
    process_pdb_files(data_file_paths, args.output, args.processes,
                      args.verbose)

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
    engine = create_engine(f'sqlite:///{output_path}', echo=False)
    BASE.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    total_pdbs = 0
    succeeded = 0
    failed = 0
    total_loops = 0
    n = 100
    data_chunks = [data_file_paths[i:i + n]
                   for i in range(0, len(data_file_paths), n)]
    for paths in data_chunks:
        if processes > 1:
            with multiprocessing.Pool(processes) as pool:
                loop_lists = pool.map(fault_tolerant_gather, paths)
        else:
            loop_lists = map(fault_tolerant_gather, paths)
        total_pdbs += len(paths)
        for loop in itertools.chain(*loop_lists):
            total_loops += 1
            if isinstance(loop, Exception):
                e = loop
                failed += 1
                if verbose:
                    print(f'Failed to process {e.path}:\n {str(e)}\n')
                else:
                    print(f'Failed to process {e.path}')
                continue
            loop_entry = Loops(**loop)
            session.add(loop_entry)
            succeeded += 1
        session.commit()
    print(f'Finished processing {total_pdbs} PDB files:\n'
          f'{succeeded} succeeded, {failed} failed, {total_loops} added to '
          f'database.')
    return


def fault_tolerant_gather(path: str) -> Union[dict, Exception]:
    """Returns the exception rather than raising on failure."""
    try:
        loops = gather_loops_from_pdb(path)
    except Exception as e:
        e.path = path
        loops = [e]
    return loops


def get_args():
    """Loads arguments for the main."""
    description = (
        "Commandline tool for creating the loop database for ISAMBARD "
        "loop modelling functionality."
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        'data_dir',
        metavar='DATA_DIRECTORY',
        help="Folder containing PDB files used to create the loop database."
    )
    parser.add_argument(
        '-o', '--output',
        help=("Path and file name for the database that's created. "
              "If no path is supplied, a file called loops.db will be "
              "created in the current working directory."),
        metavar='OUTPUT_DIRECTORY',
        default="loops.db",
        type=str
    )
    parser.add_argument(
        '-p', '--processes',
        help="Number of processes used to build the loop database.",
        type=int,
        default=1
    )
    parser.add_argument(
        '-r', '--recursive',
        help="Will recursively search the DATA_DIRECTORY for PDB files.",
        action="store_true"
    )
    parser.add_argument(
        '-v', '--verbose',
        help="Verbose output.",
        action="store_true"
    )
    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    main()
