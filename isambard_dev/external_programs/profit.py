import os
import re
import subprocess
import tempfile

from settings import global_settings


def run_profit(pdb1, pdb2, path1=False, path2=False, path_to_cmd_file=None,return_pdb_string=False,align_type=None):

    """Takes 2 PDB strings, carries out rmsd superposition using ProFit, optionally returns fitted structure as a string

    Parameters
    ----------
    pdb1 : str
        PDB as string or path
    pdb2 : str
        PDB as string or path
    path1 : bool
        Whether pdb1 is a string or filepath
    path2 : bool
        Whether pdb2 is a string or filepath
    path_to_cmd_file : None
        Optional custom command file for ProFit. Do not use if you want to use return_pdb_string=True
    return_pdb_string : bool
        Returns fitted pdb structure as a string
    align_type : None
        Used in conjunction with return_pdb_string=True and should be set to one of 'all', 'bb', or 'ca to specify
        alignment of all atoms, backbone atoms or just c-alpha atoms respectively

    Returns
    -------
    rmsds : []
        list of ca, bb and all-atom rmsds from superposition
    output_pdb : str
        (Optional) PDB string of overlaid, fitted structure (i.e., pdb2 superposed onto pdb1)

    """

    alignments = {'all': '*', 'bb': 'n,ca,c,o', 'ca': 'ca'}
    output_pdb = None
    output_file_path = None

    if path_to_cmd_file != None and return_pdb_string:
        raise ValueError("Cannot run ProFit with a custom command file and output a PDB string at the same time")

    try:
        if not path1:
            # if statement added to be sure that encode is only called on string type.
            if type(pdb1) == str:
                pdb1 = pdb1.encode()
            pdb1_tmp = tempfile.NamedTemporaryFile(delete=False)
            pdb1_tmp.write(pdb1)
            pdb1_tmp.seek(0)
            pdb1 = pdb1_tmp.name
        if not path2:
            if type(pdb2) == str:
                pdb2 = pdb2.encode()
            pdb2_tmp = tempfile.NamedTemporaryFile(delete=False)
            pdb2_tmp.write(pdb2)
            pdb2_tmp.seek(0)
            pdb2 = pdb2_tmp.name

        if path_to_cmd_file:
            cmd_file_path = path_to_cmd_file

        elif return_pdb_string:

            cmd_list = ['ignoremissing', 'align']

            if not align_type:

                cmd_list.append('atoms *')

            else:

                if align_type in alignments:

                    cmd_list.append("atoms {}".format(alignments[align_type]))
                else:
                    raise ValueError("align_type should be one of 'ca','bb' or 'all'")

            cmd_list.append('fit')

            output_file_path = tempfile.NamedTemporaryFile(delete=False)

            cmd_list.append("write {}".format(output_file_path.name))
            cmd_list.append("quit")

            tmp_cmd_file = tempfile.NamedTemporaryFile(delete=False)
            tmp_cmd_file.write(("\n".join(cmd_list)).encode())
            tmp_cmd_file.seek(0)

            cmd_file_path = tmp_cmd_file.name


        else:
            cmd_file_path = os.path.join(global_settings['package_path'], 'external_programs',
                                         'profit_cmd_files', 'all_atom_cmds.txt')
        profit_out = subprocess.check_output([global_settings['profit']['path'], '-f', cmd_file_path, pdb1, pdb2])
        rmsd_strs = re.findall('RMS: ([0-9.]+)', profit_out.decode())
        if len(rmsd_strs) != 3 and not return_pdb_string:
            raise ValueError(
                    'ProFit did not return an RMS value, check command file. See ProFit output:\n\n{}\n'
                    'PROFIT FAILED TO RUN: SEE LOG ABOVE'.format(profit_out))
        # RMSDs should contain the CA, backbone and all atom scores

        if return_pdb_string and output_file_path != None:
            output_pdb = output_file_path.read().decode()
        rmsds = [float(x) for x in rmsd_strs]

    finally:
        if not path1:
            pdb1_tmp.close()
            os.remove(pdb1_tmp.name)
        if not path2:
            pdb2_tmp.close()
            os.remove(pdb2_tmp.name)

        if return_pdb_string and output_file_path != None:
            output_file_path.close()
            tmp_cmd_file.close()
            os.remove(output_file_path.name)
            os.remove(tmp_cmd_file.name)

    if return_pdb_string and output_pdb != None:

        return rmsds, output_pdb

    else:
        return rmsds
