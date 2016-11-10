import subprocess
import tempfile
from pathlib import Path

from settings import global_settings


def run_reduce(input_file, path=True):
    """ Runs reduce on a pdb or mmol file at the specified path.

    Notes
    -----
    Runs Reduce programme to add missing protons to a PDB file.
    Reduce published in:
    Word, et al.(1999) "Asparagine and glutamine: using hydrogen atom contacts in the choice of sidechain amide
    orientation" J. Mol. Biol. 285, 1735-1747.
    Requires the reduce executable and reduce_wwPDB_het_dict.txt located in a directory specified in global_settings.
    These can be downloaded from:
    http://kinemage.biochem.duke.edu/software/reduce.php

    Parameters
    ----------
    input_file : str
        Path to file to add protons to or structure in mmol/pdb format.
    path : bool
        True if input_file is a path.

    Returns
    -------
    reduce_mmol : str
        Structure file with protons added.
    reduce_message : str
        Messages generated while running Reduce.
    """
    if path:
        input_path = Path(input_file)
        if not input_path.exists():
            print('No file found at', path)
            return None, None
    else:
        pathf = tempfile.NamedTemporaryFile()
        encoded_input = input_file.encode()
        pathf.write(encoded_input)
        pathf.seek(0)
        file_path = pathf.name
        input_path = Path(file_path)
    reduce_folder = Path(global_settings['reduce']['folder'])
    reduce_exe = reduce_folder / global_settings['reduce']['path']
    reduce_dict = reduce_folder / 'reduce_wwPDB_het_dict.txt'
    try:
        reduce_output = subprocess.run([str(reduce_exe), '-build', '-DB', str(reduce_dict), str(input_path)],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e, '\nThe Reduce executable cannot be found. Ensure the location and filename are specified in settings.')
        return None, None
    try:
        reduced_mmol = reduce_output.stdout.decode()
    except UnicodeDecodeError:
        print("Reduce could not detect any missing protons in the protein. Using the original structure.")
        if path:
            reduced_mmol = input_path.read_text()
        else:
            reduced_mmol = input_file
    reduce_message = reduce_output.stderr.decode()
    if 'could not open' in reduce_message:
        print('Caution: the Reduce connectivity dictionary could not be found. Some protons may be missing. See notes.')
    return reduced_mmol, reduce_message


def reduce_output_path(path=None, pdb_name=None):
    """ Defines where output files from Reduce will be relative to input files."""
    if not path:
        if not pdb_name:
            raise NameError("Cannot save an output for a temporary file without a PDB code specified")
        pdb_name = pdb_name.lower()
        output_path = Path(global_settings['structural_database']['path'], pdb_name[1:3].lower(), pdb_name[:4].lower(),
                           'reduce', pdb_name + '_reduced.mmol')
    else:
        input_path = Path(path)
        if len(input_path.parents) > 1:
            output_path = input_path.parents[1] / 'reduce' / (input_path.stem + '_reduced' + input_path.suffix)
        else:
            output_path = input_path.parent / (input_path.stem + '_reduced' + input_path.suffix)
    return output_path


def output_reduce(input_file, path=True, pdb_name=None, force=False):
    """ Runs Reduce on a pdb or mmol file and creates a new file with the output.

    Parameters
    ----------
    input_file : str or pathlib.Path
       Path to file to run Reduce on.
    path : bool
        True if input_file is a path.
    pdb_name : str
        PDB ID of protein. Required if providing string not path.
    force : bool
        True if existing reduce outputs should be overwritten.

    Returns
    -------
    output_path : pathlib.Path
        Location of output file.
    """
    if path:
        output_path = reduce_output_path(path=input_file)
    else:
        output_path = reduce_output_path(pdb_name=pdb_name)
    if output_path.exists() and not force:
        return output_path
    reduce_mmol, reduce_message = run_reduce(input_file, path=path)
    if not reduce_mmol:
        return None
    output_path.parent.mkdir(exist_ok=True)
    output_path.write_text(reduce_mmol)
    return output_path


def output_reduce_list(path_list, force=False):
    """ Generates structure files with protons from a list of structure files."""
    output_paths = []
    for path in path_list:
        output_path = output_reduce(path, force=force)
        if output_path:
            output_paths.append(output_path)
    return output_paths


def assembly_plus_protons(input_file, path=True, pdb_name=None, save_output=False, force_save=False):
    """ Returns an Assembly with protons added by Reduce from an input structure file.

    Notes
    -----
    Looks for a pre-existing Reduce output in the standard location before running Reduce.
    If the protein contains oligosaccharides or glycans, use reduce_correct_carbohydrates.

    Parameters
    ----------
    input_file : str or pathlib.Path
        Location of file to be converted to Assembly or PDB file as string.
    path : bool
        Whether we are looking at a file or a pdb string. Defaults to file.
    pdb_name : str
        PDB ID of protein. Required if providing string not path.
    save_output : bool
        If True will save the generated assembly.
    force_save : bool
        If True will overwrite existing reduced assembly.

    Returns
    -------
    reduced_assembly : AMPAL Assembly
        Assembly of protein with protons added by Reduce.
    """
    from ampal.pdb_parser import convert_pdb_to_ampal

    if path:
        input_path = Path(input_file)
        if not pdb_name:
            pdb_name = input_path.stem[:4]
        reduced_path = reduce_output_path(path=input_path)
        if reduced_path.exists() and not save_output and not force_save:
            reduced_assembly = convert_pdb_to_ampal(str(reduced_path), pdb_id=pdb_name)
            return reduced_assembly
    if save_output:
        reduced_path = output_reduce(input_file, path=path, pdb_name=pdb_name, force=force_save)
        reduced_assembly = convert_pdb_to_ampal(str(reduced_path), path=True)
    else:
        reduce_mmol, reduce_message = run_reduce(input_file, path=path)
        if not reduce_mmol:
            return None
        reduced_assembly = convert_pdb_to_ampal(reduce_mmol, path=False, pdb_id=pdb_name)
    return reduced_assembly


__author__ = 'Kieran L. Hudson, Gail J. Bartlett'
