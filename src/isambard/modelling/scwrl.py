"""This module provides an interface to the program Scwrl4.

The Scwrl executable must be on your path. Run the `test_scwrl`
function to determine if it is available.

For more information on Scwrl see [1].

References
----------
.. [1] Krivov GG, Shapovalov MV, and Dunbrack Jr RL (2009) "Improved
   prediction of protein side-chain conformations with SCWRL4.",
   Proteins.
"""

import os
import subprocess
import tempfile
import re

import ampal


def scwrl_available():
    """True if Scwrl is available."""
    available = False
    try:
        subprocess.check_output(['Scwrl4'], stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        available = True
    except FileNotFoundError:
        print("Scwrl4 has not been found on your path. If you have already "
              "installed Scwrl but are unsure how to add it to your path, "
              "check out this: https://stackoverflow.com/a/14638025")
    return available


def run_scwrl(pdb, sequence,
              path=True, rigid_rotamer_model=True, hydrogens=False):
    """Runs SCWRL on input PDB strong or path to PDB and a sequence string.

    Parameters
    ----------
    pdb : str
        PDB string or a path to a PDB file.
    sequence : str
        Amino acid sequence for SCWRL to pack in single-letter code.
    path : bool, optional
        True if pdb is a path.
    rigid_rotamer_model : bool, optional
        If True, Scwrl will use the rigid-rotamer model, which is
        faster but less accurate.
    hydrogens : bool, optional
        If False, the hydrogens produced by Scwrl will be ommitted.

    Returns
    -------
    scwrl_std_out : str
        Std out from SCWRL.
    scwrl_pdb : str
        String of packed SCWRL PDB.

    Raises
    ------
    ChildProcessError
        Raised if SCWRL failed to run.
    """
    if path:
        with open(pdb, 'r') as inf:
            pdb = inf.read()
    pdb = pdb.encode()
    sequence = sequence.encode()
    try:
        with tempfile.NamedTemporaryFile(delete=False) as scwrl_tmp,\
                tempfile.NamedTemporaryFile(delete=False) as scwrl_seq,\
                tempfile.NamedTemporaryFile(delete=False) as scwrl_out:
            scwrl_tmp.write(pdb)
            scwrl_tmp.seek(0)  # Resets the buffer back to the first line
            scwrl_seq.write(sequence)
            scwrl_seq.seek(0)
            scwrl_command = ['Scwrl4',
                             '-i', scwrl_tmp.name,
                             '-o', scwrl_out.name,
                             '-s', scwrl_seq.name]
            if rigid_rotamer_model:
                scwrl_command.append('-v')
            if not hydrogens:
                scwrl_command.append('-h')
            scwrl_std_out = subprocess.check_output(scwrl_command)
            scwrl_out.seek(0)
            scwrl_pdb = scwrl_out.read()
    finally:
        os.remove(scwrl_tmp.name)
        os.remove(scwrl_out.name)
        os.remove(scwrl_seq.name)
    if not scwrl_pdb:
        raise ChildProcessError(
            'SCWRL failed to run. SCWRL:\n{}'.format(scwrl_std_out))
    return scwrl_std_out.decode(), scwrl_pdb.decode()


def parse_scwrl_out(scwrl_std_out, scwrl_pdb):
    """Parses SCWRL output and returns PDB and SCWRL score.
    Parameters
    ----------
    scwrl_std_out : str
        Std out from SCWRL.
    scwrl_pdb : str
        String of packed SCWRL PDB.
    Returns
    -------
    fixed_scwrl_str : str
        String of packed SCWRL PDB, with correct PDB format.
    score : float
        SCWRL Score
    """
    score = re.findall(
        r'Total minimal energy of the graph = ([-0-9.]+)', scwrl_std_out)[0]
    # Add temperature factors to SCWRL out
    split_scwrl = scwrl_pdb.splitlines()
    fixed_scwrl = []
    for line in split_scwrl:
        if len(line) < 80:
            line += ' ' * (80 - len(line))
        if re.search(r'H?E?T?ATO?M\s+\d+.+', line):
            front = line[:61]
            temp_factor = ' 0.00'
            back = line[66:]
            fixed_scwrl.append(''.join([front, temp_factor, back]))
        else:
            fixed_scwrl.append(line)
    fixed_scwrl_str = '\n'.join(fixed_scwrl) + '\n'
    return fixed_scwrl_str, float(score)


def pack_side_chains_scwrl(assembly, sequences,
                           rigid_rotamer_model=True, hydrogens=False):
    """Packs side chains onto a protein structure.

    Parameters
    ----------
    assembly : AMPAL Assembly
        AMPAL object containing some protein structure.
    sequence : [str]
        A list of amino acid sequences in single-letter code for Scwrl to pack.
    rigid_rotamer_model : bool, optional
        If True, Scwrl will use the rigid-rotamer model, which is
        faster but less accurate.
    hydrogens : bool, optional
        If False, the hydrogens produced by Scwrl will be ommitted.

    Returns
    -------
    packed_structure : AMPAL Assembly
        A new AMPAL Assembly containing the packed structure, with
        the Scwrl score in the tags.
    """
    if not scwrl_available():
        raise ValueError('Scwrl4 is unavailable on your system path.')
    protein = [x for x in assembly if isinstance(x, ampal.Polypeptide)]
    total_seq_len = sum([len(x) for x in sequences])
    total_aa_len = sum([len(x) for x in protein])
    if total_seq_len != total_aa_len:
        raise ValueError('Total sequence length ({}) does not match '
                         'total Polypeptide length ({}).'.format(
                             total_seq_len, total_aa_len))
    if len(protein) != len(sequences):
        raise ValueError('Number of sequences ({}) does not match '
                         'number of Polypeptides ({}).'.format(
                             len(sequences), len(protein)))
    scwrl_std_out, scwrl_pdb = run_scwrl(
        assembly.pdb, ''.join(sequences), path=False,
        rigid_rotamer_model=rigid_rotamer_model, hydrogens=hydrogens)
    packed_structure, scwrl_score = parse_scwrl_out(scwrl_std_out, scwrl_pdb)
    new_assembly = ampal.load_pdb(packed_structure, path=False)
    new_assembly.tags['scwrl_score'] = scwrl_score
    return new_assembly


__author__ = 'Christopher W. Wood'
