import os
import subprocess
import tempfile
import re

from settings import global_settings
from tools.isambard_warnings import check_availability


def test_scwrl():
    is_scwrl_available = False
    if os.path.isfile(global_settings['scwrl']['path']):
        try:
            subprocess.check_output([global_settings['scwrl']['path']], stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            is_scwrl_available = True
    return is_scwrl_available


def run_scwrl(pdb, sequence, path=True):
    """Runs SCWRL on input PDB strong or path to PDB and a sequence string.

    Parameters
    ----------
    pdb : str
        PDB string or a path to a PDB file.
    sequence : str
        Amino acid sequence for SCWRL to pack in single-letter code.
    path : bool
        True if pdb is a path.

    Returns
    -------
    scwrl_std_out : str
        Std out from SCWRL.
    scwrl_pdb : str
        String of packed SCWRL PDB.

    Raises
    ------
    IOError
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
            if not global_settings['scwrl']['rigid_rotamer_model']:
                scwrl_std_out = subprocess.check_output([global_settings['scwrl']['path'],
                                                         '-i', scwrl_tmp.name,
                                                         '-o', scwrl_out.name,
                                                         '-s', scwrl_seq.name])
            else:
                scwrl_std_out = subprocess.check_output([global_settings['scwrl']['path'],
                                                         '-v',  # Rigid rotamer model
                                                         '-i', scwrl_tmp.name,
                                                         '-o', scwrl_out.name,
                                                         '-s', scwrl_seq.name])
            scwrl_out.seek(0)
            scwrl_pdb = scwrl_out.read()
    finally:
        os.remove(scwrl_tmp.name)
        os.remove(scwrl_out.name)
        os.remove(scwrl_seq.name)
    if not scwrl_pdb:
        raise IOError('SCWRL failed to run. SCWRL:\n{}'.format(scwrl_std_out))
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
    score = re.findall(r'Total minimal energy of the graph = ([-0-9.]+)', scwrl_std_out)[0]
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


@check_availability('scwrl', test_scwrl, global_settings)
def pack_sidechains(pdb, sequence, path=False):
    """Packs sidechains onto a given PDB file or string."""
    scwrl_std_out, scwrl_pdb = run_scwrl(pdb, sequence, path=path)
    return parse_scwrl_out(scwrl_std_out, scwrl_pdb)


__author__ = 'Christopher W. Wood'
