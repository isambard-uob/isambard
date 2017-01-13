import os
import subprocess
import tempfile

from settings import global_settings
from tools.isambard_warnings import check_availability


def test_dssp():
    is_dssp_available = False
    if os.path.isfile(global_settings['dssp']['path']):
        try:
            subprocess.check_output([global_settings['dssp']['path'], '--version'], stderr=subprocess.DEVNULL)
            is_dssp_available = True
        except:
            pass
    return is_dssp_available


# TODO: make the format closer to SCWRL and BUDE
# TODO: Update to APM format to remove split_pdb_lines
@check_availability('dssp', test_dssp, global_settings)
def run_dssp(pdb, path=True, outfile=None):
    """Uses DSSP to find helices and extracts helices from a pdb file or string.

    Parameters
    ----------
    pdb : str
        Path to pdb file or string.
    path : bool
        Indicates if pdb is a path or a string.
    outfile : str, optional
        Filepath for storing the dssp output.

    Returns
    -------
    dssp_out : str
        Std out from DSSP.
    """
    if not path:
        # if statement added to be sure that encode is only called on string type.
        if type(pdb) == str:
            pdb = pdb.encode()
        try:
            temp_pdb = tempfile.NamedTemporaryFile(delete=False)
            temp_pdb.write(pdb)
            temp_pdb.seek(0)
            dssp_out = subprocess.check_output([global_settings['dssp']['path'], temp_pdb.name])
            temp_pdb.close()
        finally:
            os.remove(temp_pdb.name)
    else:
        dssp_out = subprocess.check_output([global_settings['dssp']['path'], pdb])

    # Python 3 string formatting.
    dssp_out = dssp_out.decode()

    if outfile:
        with open(outfile, 'w') as outf:
            outf.write(dssp_out)
    return dssp_out


def extract_all_ss_dssp(in_dssp, path=True):
    """Uses DSSP to extract secondary structure information on every residue in the input dssp file.

    Parameters
    ----------
    in_dssp : str
        Path to DSSP file.

    Returns
    -------
    dssp_residues : [list]
        Each internal list contains:
            [0] int Residue number
            [1] str Secondary structure type
            [2] str Chain identifier
            [3] str Residue type
            [4] float Phi torsion angle
            [5] float Psi torsion angle
    """

    if path:
        with open(in_dssp, 'r') as inf:
            dssp_out = inf.read()
    else:
        dssp_out = in_dssp[:]
    dssp_residues = []
    go = False
    for line in dssp_out.splitlines():
        if go:
            try:
                res_num = int(line[5:10].strip())
                chain = line[10:12].strip()
                residue = line[13]
                ss_type = line[16]
                phi = float(line[103:109].strip())
                psi = float(line[109:116].strip())
                dssp_residues.append([res_num, ss_type, chain,
                                      residue, phi, psi])
            except ValueError:
                pass
        else:
            if line[2] == '#':
                go = True
            pass
    return dssp_residues


def extract_solvent_accessibility_dssp(in_dssp, path=True):
    """Uses DSSP to extract solvent accessibilty information on every residue in the input dssp file.

    Notes
    -----
    For more information on the solvent accessibility metrics used in dssp, see:
    http://swift.cmbi.ru.nl/gv/dssp/HTML/descrip.html#ACC
    In the dssp files value is labeled 'ACC'.

    Parameters
    ----------
    in_dssp : str
        Path to DSSP file.
    path : bool
        Indicates if in_dssp is a path or a string.

    Returns
    -------
    dssp_residues : list
        Each internal list contains:
            [0] int Residue number
            [1] str Chain identifier
            [2] str Residue type
            [3] int dssp solvent accessibilty
    """
    if path:
        with open(in_dssp, 'r') as inf:
            dssp_out = inf.read()
    else:
        dssp_out = in_dssp[:]
    dssp_residues = []
    go = False
    for line in dssp_out.splitlines():
        if go:
            try:
                res_num = int(line[5:10].strip())
                chain = line[10:12].strip()
                residue = line[13]
                acc = int(line[35:38].strip())
                # It is IMPORTANT that acc remains the final value of the returned list, due to its usage in
                # isambard.ampal.base_ampal.tag_dssp_solvent_accessibility
                dssp_residues.append([res_num, chain, residue, acc])
            except ValueError:
                pass
        else:
            if line[2] == '#':
                go = True
            pass
    return dssp_residues


def extract_helices_dssp(in_pdb):
    """Uses DSSP to find helices and extracts helices from a pdb file.

    Returns a length 3 list with a helix id, the chain id and a dict containing the coordinates of each residues CA.
    """
    from ampal.pdb_parser import split_pdb_lines

    dssp_out = subprocess.check_output([global_settings['dssp']['path'], in_pdb])
    helix = 0
    helices = []
    h_on = False
    for line in dssp_out.splitlines():
        dssp_line = line.split()
        try:
            if dssp_line[4] == 'H':
                if helix not in [x[0] for x in helices]:
                    helices.append([helix, dssp_line[2], {int(dssp_line[1]): None}])
                else:
                    helices[helix][2][int(dssp_line[1])] = None
                h_on = True
            else:
                if h_on:
                    helix += 1
                    h_on = False
        except IndexError:
            pass
    with open(in_pdb, 'r') as pdb:
        pdb_atoms = split_pdb_lines(pdb.read())
    for atom in pdb_atoms:
        for helix in helices:
            if (atom[2] == "CA") and (atom[5] == helix[1]) and (atom[6] in helix[2].keys()):
                helix[2][atom[6]] = tuple(atom[8:11])
    return helices


def extract_pp_helices(in_pdb):
    phi = -75.0
    phi_d = 29.0
    psi = 145.0
    psi_d = 29.0

    pph_dssp = subprocess.check_output([global_settings['dssp']['path'], in_pdb])
    dssp_residues = []
    go = False
    for line in pph_dssp.splitlines():
        if go:
            res_num = int(line[:5].strip())
            chain = line[11:13].strip()
            ss_type = line[16]
            phi = float(line[103:109].strip())
            psi = float(line[109:116].strip())
            dssp_residues.append((res_num, ss_type, chain, phi, psi))
        else:
            if line[2] == '#':
                go = True
            pass
    pp_chains = []
    chain = []
    ch_on = False
    for item in dssp_residues:
        if (item[1] == ' ') and (phi - phi_d < item[3] < phi + phi_d):
            chain.append(item)
            ch_on = True
        else:
            if ch_on:
                pp_chains.append(chain)
                chain = []
                ch_on = False
    pp_chains = [x for x in pp_chains if len(x) > 1]
    pp_helices = []
    with open(in_pdb, 'r') as pdb:
        pdb_atoms = split_pdb_lines(pdb.read())
    for pp_helix in pp_chains:
        chain_id = pp_helix[0][2]
        res_range = [x[0] for x in pp_helix]
        helix = []
        for atom in pdb_atoms:
            if (atom[2] == "CA") and (atom[5] == chain_id) and (atom[6] in res_range):
                helix.append(tuple(atom[8:11]))
        pp_helices.append(helix)
    return pp_helices


def find_ss_regions(dssp_residues):
    """Separates parsed DSSP data into groups based on runs of secondary structure.

    Notes
    -----
    Example: all residues in a single helix/loop/strand will be gathered into a list, then the next secondary structure
    element will be gathered into a separate list, and so on.

    Parameters
    ----------
    dssp_residues : [list]
        Each internal list contains:
            [0] int Residue number
            [1] str Secondary structure type
            [2] str Chain identifier
            [3] str Residue type
            [4] float Phi torsion angle
            [5] float Psi torsion angle

    Returns
    -------
    fragments : [[list]]
        Lists grouped in continuous regions of secondary structure. Innermost list has the same format as above.
    """

    loops = [' ', 'B', 'S', 'T']
    current_ele = None
    fragment = []
    fragments = []
    first = True
    for ele in dssp_residues:
        if first:
            first = False
            fragment.append(ele)
        elif current_ele in loops:
            if ele[1] in loops:
                fragment.append(ele)
            else:
                fragments.append(fragment)
                fragment = [ele]
        else:
            if ele[1] == current_ele:
                fragment.append(ele)
            else:
                fragments.append(fragment)
                fragment = [ele]
        current_ele = ele[1]
    return fragments
