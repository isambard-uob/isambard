import json
import os
import re
import shutil
from pathlib import Path

import requests
import tempfile
import gzip

from ampal import AmpalContainer, Assembly, Polypeptide
from ampal.pdb_parser import convert_pdb_to_ampal
from external_programs.dssp import run_dssp
from settings import global_settings

_mmols_numbers_path = os.path.join(global_settings['package_path'], 'add_ons', 'mmols.json')
if os.path.exists(_mmols_numbers_path):
    with open(_mmols_numbers_path, 'r') as inf:
        mmols_numbers = json.loads(inf.read())
else:
    mmols_numbers = None


class FileSystem:
    """Get file paths for all files relevant to a given PDB code. Generate files if not already present.

    Check files in turn. If they are not yet present in the database folder, download or generate them as necessary.
    There are multiple mmol/dssp files generated.

    Parameters
    ----------
    code : str
        A PDB code
    data_dir : str
        Main folder within which FileSystem attributes will be stored.
        Default is the database_dir referenced in settings.json

    Attributes
    ----------
    code : str
        PDB code
    parent_dir : str
        Filepath for folder in which all structural subfolders will be stored.
    number_of_mmols : int
        Number of .mmol files associated with code in the PDBE.
    preferred_mmol : int
        The mmol number of preferred biological assembly as listed in the PDBe.
    """
    def __init__(self, code, data_dir=None):
        self.code = code
        if not data_dir:
            data_dir = global_settings["structural_database"]["path"]
        self._data_dir = data_dir
        self.parent_dir = os.path.join(self._data_dir, self.code[1:3], self.code)
        try:
            self.number_of_mmols = number_of_mmols(self.code)
        except ValueError:
            raise ValueError('No mmols files available in PDBE for {0}.'.format(self.code))
        self.preferred_mmol = preferred_mmol(self.code)

    def __repr__(self):
        return '<FileSystem for folder {0}/>'.format(self.parent_dir)

    @property
    def mmols(self):
        """ Dict of filepaths for all mmol files associated with code.

        Notes
        -----
        Downloads mmol files if not already present.

        Returns
        -------
        mmols_dict : dict, or None.
            Keys : int
                mmol number
            Values : str
                Filepath for the corresponding mmol file.
        """
        mmols_dict = {}
        mmol_dir = os.path.join(self.parent_dir, 'structures')
        if not os.path.exists(mmol_dir):
            os.makedirs(mmol_dir)
        mmol_file_names = ['{0}_{1}.mmol'.format(self.code, i) for i in range(1, self.number_of_mmols + 1)]
        mmol_files = [os.path.join(mmol_dir, x) for x in mmol_file_names]
        for i, mmol_file in enumerate(mmol_files):
            mmols_dict[i + 1] = mmol_file
            # If file does not exist yet, download the mmol and write to mmol_file.
            if not os.path.exists(mmol_file):
                get_mmol(self.code, mmol_number=i + 1, outfile=mmol_file)
        return mmols_dict

    @property
    def dssps(self):
        """ Dict of filepaths for all dssp files associated with code.

        Notes
        -----
        Runs dssp and stores writes output to files if not already present.
        Also downloads mmol files if not already present.
        Calls isambard.external_programs.dssp and so needs dssp to be installed.

        Returns
        -------
        dssps_dict : dict, or None.
            Keys : int
                mmol number
            Values : str
                Filepath for the corresponding dssp file.

        Raises
        ------
        Warning
            If any of the dssp files are empty.
        """
        dssps_dict = {}
        dssp_dir = os.path.join(self.parent_dir, 'dssp')
        if not os.path.exists(dssp_dir):
            os.makedirs(dssp_dir)
        for i, mmol_file in self.mmols.items():
            dssp_file_name = '{0}.dssp'.format(os.path.basename(mmol_file))
            dssp_file = os.path.join(dssp_dir, dssp_file_name)
            if not os.path.exists(dssp_file):
                dssp_out = run_dssp(pdb=mmol_file, path=True, outfile=dssp_file)
                if len(dssp_out) == 0:
                    raise Warning("dssp file {0} is empty".format(dssp_file))
            dssps_dict[i] = dssp_file
        return dssps_dict

    @property
    def fastas(self, download=False):
        """ Dict of filepaths for all fasta files associated with code.

        Parameters
        ----------
        download : bool
            If True, downloads the fasta file from the PDB.
            If False, uses the ampal Protein.fasta property
            Defaults to False - this is definitely the recommended behaviour.

        Notes
        -----
        Calls self.mmols, and so downloads mmol files if not already present.
        See .fasta property of isambard.ampal.base_ampal.Protein for more information.

        Returns
        -------
        fastas_dict : dict, or None.
            Keys : int
                mmol number
            Values : str
                Filepath for the corresponding fasta file.
        """
        fastas_dict = {}
        fasta_dir = os.path.join(self.parent_dir, 'fasta')
        if not os.path.exists(fasta_dir):
            os.makedirs(fasta_dir)
        for i, mmol_file in self.mmols.items():
            mmol_name = os.path.basename(mmol_file)
            fasta_file_name = '{0}.fasta'.format(mmol_name)
            fasta_file = os.path.join(fasta_dir, fasta_file_name)
            if not os.path.exists(fasta_file):
                if download:
                    pdb_url = "http://www.rcsb.org/pdb/files/fasta.txt?structureIdList={0}".format(self.code.upper())
                    r = requests.get(pdb_url)
                    if r.status_code == 200:
                        fasta_string = r.text
                    else:
                        fasta_string = None
                else:
                    a = convert_pdb_to_ampal(mmol_file)
                    # take first object if AmpalContainer (i.e. NMR structure).
                    if type(a) == AmpalContainer:
                        a = a[0]
                    fasta_string = a.fasta
                with open(fasta_file, 'w') as foo:
                    foo.write(fasta_string)
            fastas_dict[i] = fasta_file
        return fastas_dict

    @property
    def cifs(self):
        cif_dict = {}
        cif_dir = os.path.join(self.parent_dir, 'structures')
        if not os.path.exists(cif_dir):
            os.makedirs(cif_dir)
        for i in range(self.number_of_mmols):
            mmol_number = i + 1
            mmcif_file_name = '{0}_{1}.cif'.format(self.code, mmol_number)
            mmcif_file = os.path.join(cif_dir, mmcif_file_name)
            if not os.path.exists(mmcif_file):
                get_cif(code=self.code, mmol_number=mmol_number, outfile=mmcif_file)
            cif_dict[mmol_number] = mmcif_file
        return cif_dict

    @property
    def mmcif(self):
        """ Filepath for mmcif file associated with code.

        Notes
        -----
        Downloads mmcif file if not already present.

        Returns
        -------
        mmcif_file : str
            Filepath for the mmcif file.
        """
        mmcif_dir = os.path.join(self.parent_dir, 'mmcif')
        if not os.path.exists(mmcif_dir):
            os.makedirs(mmcif_dir)
        mmcif_file_name = '{0}.cif'.format(self.code)
        mmcif_file = os.path.join(mmcif_dir, mmcif_file_name)
        if not os.path.exists(mmcif_file):
            get_mmcif(code=self.code, outfile=mmcif_file)
        return mmcif_file


def is_obsolete(code):
    url = 'http://www.ebi.ac.uk/pdbe/entry/pdb/{0}'.format(code)
    r = requests.get(url)
    obs = False
    if r.status_code == 200:
        if "Entry has been obsoleted (OBS)" in r.content.decode():
            obs = True
    return obs


def number_of_mmols(code):
    """ Number of .mmol files associated with code in the PDBE.

    Notes
    -----
    This function makes a series of calls to the PDBE website using the requests module. This can make it slow!

    Parameters
    ----------
    code : str
        PDB code.

    Returns
    -------
    num_mmols : int

    Raises
    ------
    ValueError
        If no .mmol files are found at all.
        Could be due to erroneous input argument, or a problem with connecting to the PDBE.
    """
    # If num_mmols is already known, return it
    if mmols_numbers:
        if code in mmols_numbers.keys():
            mmol = mmols_numbers[code][0]
            return mmol
    counter = 1
    while True:
        pdbe_url = "http://www.ebi.ac.uk/pdbe/static/entry/download/{0}-assembly-{1}.cif.gz".format(code, counter)
        r = requests.get(pdbe_url)
        if r.status_code == 200:
            counter += 1
        else:
            break
    if counter == 1:
        while True:
            pdb_url = "http://www.rcsb.org/pdb/files/{0}.pdb{1}.gz".format(code.upper(), counter)
            r = requests.get(pdb_url)
            if r.status_code == 200 and r.encoding is None:
                counter += 1
            else:
                break
    if counter == 1:
        pdb_url = "http://files.rcsb.org/download/{0}.pdb".format(code.upper())
        r = requests.get(pdb_url)
        if r.status_code == 200:
            counter += 1
    num_mmols = counter - 1
    if num_mmols == 0:
        raise ValueError('Could not access ANY .mmol files for {0}'.format(code))
    return num_mmols


def get_mmol(code, mmol_number=None, outfile=None):
    """ Get mmol file from PDBe and return its content as a string. Write to file if outfile given.

    Parameters
    ----------
    code : str
        PDB code.
    mmol_number : int
        mmol number (biological assembly number) of file to download. Numbers from PDBe.
        If None, defaults to the preferred biological assembly listed for code on the PDBe.
    outfile : str
        Filepath. Writes returned value to this file.

    Returns
    -------
    mmol_string : str, or None
        Content of the mmol file as a string.
        None if there are no pdbe files to download, as determined by pdbe_status_code().
        None if unable to download the mmol_file from the pdbe site.

    Raises
    ------
    ValueError
        If the number of mmols for code is stored in mmols_numbers and if mmol_number is larger than this value.
    """
    if not mmol_number:
        try:
            mmol_number = preferred_mmol(code=code)
        except (ValueError, TypeError, IOError):
            print("No mmols for {0}".format(code))
            return None
    # sanity check
    if mmols_numbers:
        if code in mmols_numbers.keys():
            num_mmols = mmols_numbers[code][0]
            if mmol_number > num_mmols:
                raise ValueError('There are only {0} mmols for code {1}. mmol_number {2} is too big'
                                 .format(num_mmols, code, mmol_number))
    # Download mmol file from the PDBE webserver.
    pdbe_url = "http://www.ebi.ac.uk/pdbe/entry-files/download/{0}_{1}.mmol".format(code, mmol_number)
    r = requests.get(pdbe_url)
    if r.status_code == 200:
        mmol_string = r.text
    else:
        # Download gz pdb file from the PDB.
        pdb_url = "http://www.rcsb.org/pdb/files/{0}.pdb{1}.gz".format(code.upper(), mmol_number)
        r = requests.get(pdb_url)
        if r.status_code == 200:
            temp_gz = tempfile.NamedTemporaryFile()
            temp_gz.write(r.content)
            with gzip.open(temp_gz.name, 'rb') as foo:
                mmol_string = foo.read().decode()
        else:
            print("Could not download mmol file for {0}.\n Got requests status_code {1}".format(code, r.status_code))
            return None
    # Write to file
    if outfile and mmol_string:
        with open(outfile, 'w') as foo:
            foo.write(mmol_string)

    return mmol_string


def get_cif(code, mmol_number, outfile=None):
    """
    Parameters
    ----------
    code : str
        PDB code.
    mmol_number : int
        mmol number (biological assembly number) of file to download. Numbers from PDBe.
        If None, defaults to the preferred biological assembly listed for code on the PDBe.
    outfile : str
        Filepath. Writes returned value to this file.

    Returns
    -------
    cif_string : str, or None
        Content of the cif file as a string.
        None if unable to download the cif_file from the pdbe site.
    """
    pdbe_url = "http://www.ebi.ac.uk/pdbe/static/entry/download/{0}-assembly-{1}.cif.gz".format(code, mmol_number)
    r = requests.get(pdbe_url)
    if r.status_code == 200:
        temp_gz = tempfile.NamedTemporaryFile()
        temp_gz.write(r.content)
        with gzip.open(temp_gz.name, 'rb') as foo:
            cif_string = foo.read().decode()
    else:
        print("Could not download cif file for {0}".format(code))
        return None
    # Write to file.
    if outfile and cif_string:
        with open(outfile, 'w') as foo:
            foo.write(cif_string)
    return cif_string


def get_mmcif(code, outfile=None):
    """ Get mmcif file associated with code from PDBE.

    Parameters
    ----------
    code : str
        PDB code.
    outfile : str
        Filepath. Writes returned value to this file.

    Returns
    -------
    mmcif_file : str
        Filepath to the mmcif file.
    """
    pdbe_url = "http://www.ebi.ac.uk/pdbe/entry-files/download/{0}.cif".format(code)
    r = requests.get(pdbe_url)
    if r.status_code == 200:
        mmcif_string = r.text
    else:
        print("Could not download mmcif file for {0}".format(code))
        mmcif_string = None

    # Write to file.
    if outfile and mmcif_string:
        with open(outfile, 'w') as foo:
            foo.write(mmcif_string)

    return mmcif_string


def pdbe_status_code(code):
    """Check if a PDB code has structure files on the PDBE site.

    Parameters
    ----------
    code : str
        PDB code to check for on PDBE.

    Returns
    -------
    status_code : int
        HTTP status code of PDBE url associated with input code.
    """
    url = 'http://www.ebi.ac.uk/pdbe/entry-files/download/{0}_1.mmol'.format(code)
    r = requests.head(url=url)
    return r.status_code


def preferred_mmol(code):
    """ Get mmol number of preferred biological assembly as listed in the PDBe.

    Notes
    -----
    First checks for code in mmols.json.
    If code not yet in this json dictionary, uses requests module to scrape the PDBE for the preferred mmol number.

    Parameters
    ----------
    code : str
        A PDB code.

    Returns
    -------
    mmol : int
        mmol number of preferred assembly.

    Raises
    ------
    TypeError
        If 'mmol number' scraped is not an integer.
    """
    # If preferred mmol number is already known, return it
    if code in mmols_numbers.keys():
        mmol = mmols_numbers[code][1]
        return mmol
    elif is_obsolete(code):
        raise ValueError('Obsolete PDB code {0}'.format(code))
    # Otherwise, use requests to scrape the PDBE.
    else:
        url_string = "http://www.ebi.ac.uk/pdbe/entry/pdb/{0}/analysis".format(code)
        r = requests.get(url_string)
        if not r.ok:
            raise IOError("Could not get to url {0}".format(url_string))
        r_content = r.text
        ass = re.findall('Assembly\s\d+\s\(preferred\)', r_content)
        if len(ass) != 1:
            # To catch a strange error in the pdbe where preferred assembly is not numbered. See for example
            # http://www.ebi.ac.uk/pdbe/entry/pdb/7msi/analysis
            ass = re.findall('Assembly\s+\(preferred\)', r_content)
            if len(ass) == 1:
                return 1
            obs = re.findall('Entry has been obsoleted and replaced by another entry \(OBS\)', r_content)
            if len(obs) == 1:
                rep = re.findall('by entry <a href="/pdbe/entry/pdb/\w{4}', r_content)
                if len(rep) == 1:
                    rep = rep[0][-4:]
                    raise IOError("{0} is obsolete and has been replaced by {1}.".format(code, rep))
            raise ValueError("More than one match to preferred assembly")
        mmol = ass[0].split()[1]
        try:
            mmol = int(mmol)
        except TypeError:
            raise TypeError("Unexpected match: non-integer mmol")
        return mmol


# PDB management functions. Possibly to be moved to their own module.
def current_codes_from_pdb():
    """ Get list of all PDB codes currently listed in the PDB.

    Returns
    -------
    pdb_codes : list(str)
        List of PDB codes (in lower case).
    """
    url = 'http://www.rcsb.org/pdb/rest/getCurrent'
    r = requests.get(url)
    if r.status_code == 200:
        pdb_codes = [x.lower() for x in r.text.split('"') if len(x) == 4]
    else:
        print('Request for {0} failed with status code {1}'.format(url, r.status_code))
        return
    return pdb_codes


def obsolete_codes_from_pdb():
    """ Get list of all PDB codes listed in the PDB as being obsolete.

    Returns
    -------
    pdb_codes : list(str)
        List of PDB codes (in lower case).
    """
    url = 'http://www.rcsb.org/pdb/rest/getObsolete'
    r = requests.get(url)
    if r.status_code == 200:
        pdb_codes = [x.lower() for x in r.text.split('"') if len(x) == 4]
    else:
        print('Request for {0} failed with status code {1}'.format(url, r.status_code))
        return
    return pdb_codes


def local_pdb_codes(data_dir=None):
    """ Get list of PDB codes stored in a folder (FileSystem folder hierarchy expected within data_dir).

    If no folder is specified, use the database_dir defined in settings.json.

    Parameters
    ----------
    data_dir: str
        Filepath to a folder containing the PDB folder hierarchy (eg data_dir/eb/2ebo)

    Returns
    -------
    pdb_code_list : list(str)
        PDB codes present in data_dir.

    """
    if not data_dir:
        data_dir = global_settings["structural_database"]["path"]
    p = Path(data_dir)
    pdb_parent_dirs = [x for x in p.iterdir() if x.is_dir() and len(x.parts[-1]) == 2]
    pdb_folders = [x for test in pdb_parent_dirs for x in test.iterdir() if x.is_dir()]
    pdb_code_list = [x.parts[-1] for x in pdb_folders if len(x.parts[-1]) == 4]
    return pdb_code_list


def make_code_obsolete(code):
    """ Moves folders associated with PDB code to obsolete folder in global_settings["database_dir"]

    Parameters
    ----------
    code : str
        PDB accession code

    Returns
    -------
    None
    """
    fs = FileSystem(code=code)
    if os.path.exists(fs.parent_dir):
        # Move to obsolete folder
        destination_dir = os.path.join(fs._data_dir, 'obsolete', code[1:3], code)
        if os.path.exists(destination_dir):
            shutil.rmtree(destination_dir)
        shutil.move(fs.parent_dir, destination_dir)

        # Remove containing (two-letter) folder if empty, else pass.
        two_letter_dir = os.path.dirname(fs.parent_dir)
        try:
            os.rmdir(two_letter_dir)
        except OSError:
            pass

    return

__author__ = 'Jack W. Heal'
