import json
import requests
import os

from settings import global_settings

_amino_acids_json_path = os.path.join(global_settings['package_path'], 'tools', 'amino_acids.json')
with open(_amino_acids_json_path, 'r') as inf:
    amino_acids_dict = json.loads(inf.read())

water_mass = 18.01528

ideal_backbone_bond_lengths = {
    # Ideal bond distances from:
    #  Schulz, G. E, and R. Heiner Schirmer. Principles Of Protein Structure. New York: Springer-Verlag, 1979.
    'n_ca': 1.47,
    'ca_c': 1.53,
    'c_o': 1.24,
    # peptide bond length for residues.
    'c_n': 1.32,
}

ideal_backbone_bond_angles = {
    # Ideal bond angles from:
    #  Schulz, G. E, and R. Heiner Schirmer. Principles Of Protein Structure. New York: Springer-Verlag, 1979.
    'trans': {
        'n_ca_c': 110.0,
        'ca_c_o': 121.0,
        'ca_c_n': 114.0,
        'c_n_ca': 123.0,
        'o_c_n': 125.0,
    },
    'cis': {
        'n_ca_c': 110.0,
        'ca_c_o': 119.0,
        'ca_c_n': 118.0,
        'c_n_ca': 126.0,
        'o_c_n': 123.0
    }
}

residue_mwt = {
    'A': 71.0779, 'C': 103.1429, 'D': 115.0874, 'E': 129.114, 'F': 147.1739,
    'G': 57.0513, 'H': 137.1393, 'I': 113.1576, 'K': 128.1723, 'L': 113.1576,
    'M': 131.1961, 'N': 114.1026, 'P': 97.1152, 'Q': 128.1292, 'R': 156.1857,
    'S': 87.0773, 'T': 101.1039, 'V': 99.1311, 'W': 186.2099, 'Y': 163.1733,
    'X': 57.0513
    }

residue_charge = {
    'A': 0, 'C': -1, 'D': -1, 'E': -1, 'F': 0,
    'G': 0, 'H': +1, 'I': 0, 'K': +1, 'L': 0,
    'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': +1,
    'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': -1,
    'N-term': +1, 'C-term': -1, 'X': 0
    }

residue_pka = {
    'A': 0.0, 'C': 8.3, 'D': 3.65, 'E': 4.25, 'F': 0.0,
    'G': 0.0, 'H': 6.1, 'I': 0.0, 'K': 10.53, 'L': 0.0,
    'M': 0.0, 'N': 0.0, 'P': 0.0, 'Q': 0.0, 'R': 12.48,
    'S': 0.0, 'T': 0.0, 'V': 0.0, 'W': 0.0, 'Y': 10.1,
    'N-term': 8.0, 'C-term': 3.1, 'X': 0.0
    }

residue_ext_280 = {
    'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
    'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
    'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0,
    'S': 0, 'T': 0, 'V': 0, 'W': 5690, 'Y': 1280,
    'X': 0
    }

standard_amino_acids = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }

side_chain_dihedrals = {
        'ARG': [
            ['N', 'CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
            ['CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
            ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
            ['CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2']],
        'ASN': [
            ['N', 'CA', 'CB', 'CG', 'OD1', 'ND2'],
            ['CA', 'CB', 'CG', 'OD1', 'ND2']],
        'ASP': [
            ['N', 'CA', 'CB', 'CG', 'OD1', 'OD2'],
            ['CA', 'CB', 'CG', 'OD1', 'OD2']],
        'CYS': [['N', 'CA', 'CB', 'SG']],
        'GLN': [
            ['N', 'CA', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
            ['CA', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
            ['CB', 'CG', 'CD', 'OE1', 'NE2']],
        'GLU': [
            ['N', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
            ['CA', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
            ['CB', 'CG', 'CD', 'OE1', 'OE2']],
        'HIS': [
            ['N', 'CA', 'CB', 'CG', 'ND1', 'CE1', 'ND2'],
            ['CA', 'CB', 'CG', 'ND1', 'CE1', 'ND2']],
        'ILE': [
            ['N', 'CA', 'CB', 'CG1', 'CG2', 'CD1'],
            ['CA', 'CB', 'CG1', 'CD1', 'CG2']],
        'LEU': [
            ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2'],
            ['CA', 'CB', 'CG', 'CD1', 'CD2']],
        'LYS': [
            ['N', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ'],
            ['CA', 'CB', 'CG', 'CD', 'CE', 'NZ'],
            ['CB', 'CG', 'CD', 'CE', 'NZ'],
            ['CG', 'CD', 'CE', 'NZ']],
        'MET': [
            ['N', 'CA', 'CB', 'CG', 'SD', 'CE'],
            ['CA', 'CB', 'CG', 'SD', 'CE'],
            ['CB', 'CG', 'SD', 'CE']],
        'PHE': [
            ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']],
        'PRO': [
            ['N', 'CA', 'CB', 'CG', 'CD'],
            ['CA', 'CB', 'CG', 'CD']],
        'SER': [['N', 'CA', 'CB', 'OG']],
        'THR': [['N', 'CA', 'CB', 'OG1', 'CG2']],
        'TRP': [
            ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'NE1', 'CZ2', 'CZ3', 'CH2'],
            ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'NE1', 'CZ2', 'CZ3', 'CH2']],
        'TYR': [
            ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
            ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']],
        'VAL': [['N', 'CA', 'CB', 'CG1', 'CG2']]
    }

# Data taken from http://web.expasy.org/protscale/ unless otherwise stated. Original reference also given.
# Levitt M. Biochemistry 17:4277-4285(1978)
a_helix_Levitt = {
    'A': 1.29, 'C': 1.11, 'D': 1.04, 'E': 1.44, 'F': 1.07,
    'G': 0.56, 'H': 1.22, 'I': 0.97, 'K': 1.23, 'L': 1.3,
    'M': 1.47, 'N': 0.9, 'P': 0.52, 'Q': 1.27, 'R': 0.96,
    'S': 0.82, 'T': 0.82, 'V': 0.91, 'W': 0.99, 'Y': 0.72
    }


# Janin J. Nature 277:491-492(1979)
accessibility_Janin = {
    'A': 6.6, 'C': 0.9, 'D': 7.7, 'E': 5.7, 'F': 2.4,
    'G': 6.7, 'H': 2.5, 'I': 2.8, 'K': 10.3, 'L': 4.8,
    'M': 1.0, 'N': 6.7, 'P': 4.8, 'Q': 5.2, 'R': 4.5,
    'S': 9.4, 'T': 7.0, 'V': 4.5, 'W': 1.4, 'Y': 5.1
    }

# Bhaskaran R., Ponnuswamy P.K. Int. J. Pept. Protein. Res. 32:242-255(1988)
avg_flex_index = {
    'A': 0.36, 'C': 0.35, 'D': 0.51, 'E': 0.5, 'F': 0.31,
    'G': 0.54, 'H': 0.32, 'I': 0.46, 'K': 0.47, 'L': 0.37,
    'M': 0.3, 'N': 0.46, 'P': 0.51, 'Q': 0.49, 'R': 0.53,
    'S': 0.51, 'T': 0.44, 'V': 0.39, 'W': 0.31, 'Y': 0.42
    }

# Levitt M. Biochemistry 17:4277-4285(1978)
beta_sheet_Levitt = {
    'A': 0.9, 'C': 0.74, 'D': 0.72, 'E': 0.75, 'F': 1.32,
    'G': 0.92, 'H': 1.08, 'I': 1.45, 'K': 0.77, 'L': 1.02,
    'M': 0.97, 'N': 0.76, 'P': 0.64, 'Q': 0.8, 'R': 0.99,
    'S': 0.95, 'T': 1.21, 'V': 1.49, 'W': 1.14, 'Y': 1.25
    }

# Levitt M. Biochemistry 17:4277-4285(1978)
beta_turn_Levitt = {
    'A': 0.77, 'C': 0.81, 'D': 1.41, 'E': 0.99, 'F': 0.59,
    'G': 1.64, 'H': 0.68, 'I': 0.51, 'K': 0.96, 'L': 0.58,
    'M': 0.41, 'N': 1.28, 'P': 1.91, 'Q': 0.98, 'R': 0.88,
    'S': 1.32, 'T': 1.04, 'V': 0.47, 'W': 0.76, 'Y': 1.05
    }

# Zimmerman J.M., Eliezer N., Simha R. J. Theor. Biol. 21:170-201(1968)
bulkiness = {
    'A': 11.5, 'C': 13.46, 'D': 11.68, 'E': 13.57, 'F': 19.8,
    'G': 3.4, 'H': 13.69, 'I': 21.4, 'K': 15.71, 'L': 21.4,
    'M': 16.25, 'N': 12.82, 'P': 17.43, 'Q': 14.45, 'R': 14.28,
    'S': 9.47, 'T': 15.77, 'V': 21.57, 'W': 21.67, 'Y': 18.03
    }

# Kyte J., Doolittle R.F. J. Mol. Biol. 157:105-132(1982)
hydropathicity = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }

# http://biopython.org/DIST/docs/api/Bio.PDB.DSSP%27-pysrc.html
# Sander & Rost, (1994), Proteins, 20:216-226
max_asa = {
    'A': 106.0, 'C': 135.0, 'D': 163.0, 'E': 194.0, 'F': 197.0,
    'G': 84.0, 'H': 184.0, 'I': 169.0, 'K': 205.0, 'L': 164.0,
    'M': 188.0, 'N': 157.0, 'P': 136.0, 'Q': 198.0, 'R': 135.0,
    'S': 130.0, 'T': 142.0, 'V': 142.0, 'W': 227.0, 'Y': 222.0
    }

# Grantham R. Science 185:862-864(1974)
polarity_Grantham = {
    'A': 8.1, 'C': 5.5, 'D': 13.0, 'E': 12.3, 'F': 5.2,
    'G': 9.0, 'H': 10.4, 'I': 5.2, 'K': 11.3, 'L': 4.9,
    'M': 5.7, 'N': 11.6, 'P': 8.0, 'Q': 10.5, 'R': 10.5,
    'S': 9.2, 'T': 8.6, 'V': 5.9, 'W': 5.4, 'Y': 6.2
    }

# Zimmerman J.M., Eliezer N., Simha R. J. Theor. Biol. 21:170-201(1968)
polarity_Zimmerman = {
    'A': 0.0, 'C': 1.48, 'D': 49.7, 'E': 49.9, 'F': 0.35,
    'G': 0.0, 'H': 51.6, 'I': 0.13, 'K': 49.5, 'L': 0.13,
    'M': 1.43, 'N': 3.38, 'P': 1.58, 'Q': 3.53, 'R': 52.0,
    'S': 1.67, 'T': 1.66, 'V': 0.13, 'W': 2.1, 'Y': 1.61
    }

# Fraga S. Can. J. Chem. 60:2606-2610(1982)
recognition_factors = {
    'A': 78.0, 'C': 89.0, 'D': 81.0, 'E': 78.0, 'F': 81.0,
    'G': 84.0, 'H': 84.0, 'I': 88.0, 'K': 87.0, 'L': 85.0,
    'M': 80.0, 'N': 94.0, 'P': 91.0, 'Q': 87.0, 'R': 95.0,
    'S': 107.0, 'T': 93.0, 'V': 89.0, 'W': 104.0, 'Y': 84.0
    }

# Jones. D.D. J. Theor. Biol. 50:167-184(1975)
refractivity = {
    'A': 4.34, 'C': 35.77, 'D': 12.0, 'E': 17.26, 'F': 29.4,
    'G': 0.0, 'H': 21.81, 'I': 19.06, 'K': 21.29, 'L': 18.78,
    'M': 21.64, 'N': 13.28, 'P': 10.93, 'Q': 17.56, 'R': 26.66,
    'S': 6.35, 'T': 11.01, 'V': 13.92, 'W': 42.53, 'Y': 31.53
    }

# Dayhoff M.O., Schwartz R.M., Orcutt B.C. In "Atlas of Protein Sequence and Structure", Vol.5, Suppl.3 (1978)
relative_mutability = {
    'A': 100.0, 'C': 20.0, 'D': 106.0, 'E': 102.0, 'F': 41.0,
    'G': 49.0, 'H': 66.0, 'I': 96.0, 'K': 56.0, 'L': 40.0,
    'M': 94.0, 'N': 134.0, 'P': 56.0, 'Q': 93.0, 'R': 65.0,
    'S': 120.0, 'T': 97.0, 'V': 74.0, 'W': 18.0, 'Y': 41.0
    }

# Meek J.L. Proc. Natl. Acad. Sci. USA 77:1632-1636(1980)
retention_coeff_hplc_pH7pt4 = {
    'A': 0.5, 'C': -6.8, 'D': -8.2, 'E': -16.9, 'F': 13.2,
    'G': 0.0, 'H': -3.5, 'I': 13.9, 'K': 0.1, 'L': 8.8,
    'M': 4.8, 'N': 0.8, 'P': 6.1, 'Q': -4.8, 'R': 0.8,
    'S': 1.2, 'T': 2.7, 'V': 2.7, 'W': 14.9, 'Y': 6.1
    }

# Zhao, G., London E. Protein Sci. 15:1987-2001(2006)
transmembrane_tendancy = {
    'A': 0.38, 'C': -0.3, 'D': -3.27, 'E': -2.9, 'F': 1.98,
    'G': -0.19, 'H': -1.44, 'I': 1.97, 'K': -3.46, 'L': 1.82,
    'M': 1.4, 'N': -1.62, 'P': -1.44, 'Q': -1.84, 'R': -2.57,
    'S': -0.53, 'T': -0.32, 'V': 1.46, 'W': 1.53, 'Y': 0.49
    }

# Bairoch A. Release notes for UniProtKB/Swiss-Prot release 2013_04 - April 2013
uniprot_composition_2013 = {
    'A': 8.25, 'C': 1.37, 'D': 5.45, 'E': 6.75, 'F': 3.86,
    'G': 7.07, 'H': 2.27, 'I': 5.96, 'K': 5.84, 'L': 9.66,
    'M': 2.42, 'N': 4.06, 'P': 4.7, 'Q': 3.93, 'R': 5.53,
    'S': 6.56, 'T': 5.34, 'V': 6.87, 'W': 1.08, 'Y': 2.92
    }


number_of_codons = {
    'A': 4, 'C': 1, 'D': 2, 'E': 2, 'F': 2,
    'G': 4, 'H': 2, 'I': 3, 'K': 2, 'L': 6,
    'M': 1, 'N': 2, 'P': 4, 'Q': 2, 'R': 6,
    'S': 6, 'T': 4, 'V': 4, 'W': 1, 'Y': 2
    }

# pI, pk_COOH, pK_NH3, pK_Rgroup all taken from:
# http://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html#prop'
# D. R. Lide, Handbook of Chemistry and Physics, 72nd Edition, CRC Press, Boca Raton, FL, 1991.
pI = {
    'A': 6.0, 'C': 5.07, 'D': 2.77, 'E': 3.22, 'F': 5.48,
    'G': 5.97, 'H': 7.59, 'I': 6.02, 'K': 9.74, 'L': 5.98,
    'M': 5.74, 'N': 5.41, 'P': 6.3, 'Q': 5.65, 'R': 10.76,
    'S': 5.68, 'T': 5.6, 'V': 5.96, 'W': 5.89, 'Y': 5.66
    }

pK_COOH = {
    'A': 2.34, 'C': 1.96, 'D': 1.88, 'E': 2.19, 'F': 1.83,
    'G': 2.34, 'H': 1.82, 'I': 2.36, 'K': 2.18, 'L': 2.36,
    'M': 2.28, 'N': 2.02, 'P': 1.99, 'Q': 2.17, 'R': 2.17,
    'S': 2.21, 'T': 2.09, 'V': 2.32, 'W': 2.83, 'Y': 2.2
    }

pK_NH3 = {
    'A': 9.69, 'C': 10.28, 'D': 9.6, 'E': 9.67, 'F': 9.13,
    'G': 9.6, 'H': 9.17, 'I': 9.6, 'K': 8.95, 'L': 9.6,
    'M': 9.21, 'N': 8.8, 'P': 10.6, 'Q': 9.13, 'R': 9.04,
    'S': 9.15, 'T': 9.1, 'V': 9.62, 'W': 9.39, 'Y': 9.11
    }

pK_Rgroup = {
    'A': None, 'C': 8.18, 'D': 3.65, 'E': 4.25, 'F': None,
    'G': None, 'H': 6.0, 'I': None, 'K': 10.53, 'L': None,
    'M': None, 'N': None, 'P': None, 'Q': None, 'R': 12.48,
    'S': None, 'T': None, 'V': None, 'W': None, 'Y': 10.07
    }


def get_aa_code(aa_letter):
    """ Get three-letter aa code if possible. If not, return None.

    If three-letter code is None, will have to find this later from the filesystem.

    Parameters
    ----------
    aa_letter : str
        One-letter amino acid code.

    Returns
    -------
    aa_code : str, or None
        Three-letter aa code.
    """
    aa_code = None
    if aa_letter != 'X':
        for key, val in standard_amino_acids.items():
            if key == aa_letter:
                aa_code = val
    return aa_code


def get_aa_letter(aa_code):
    """ Get one-letter version of aa_code if possible. If not, return 'X'.

    Parameters
    ----------
    aa_code : str
        Three-letter amino acid code.

    Returns
    -------
    aa_letter : str
        One-letter aa code.
        Default value is 'X'.
    """
    aa_letter = 'X'
    for key, val in standard_amino_acids.items():
        if val == aa_code:
            aa_letter = key

    return aa_letter


def get_aa_info(code):
    """Get dictionary of information relating to a new amino acid code not currently in the database.

    Notes
    -----
    Use this function to get a dictionary that is then to be sent to the function add_amino_acid_to_json().
    use to fill in rows of amino_acid table for new amino acid code.

    Parameters
    ----------
    code : str
        Three-letter amino acid code.

    Raises
    ------
    IOError
        If unable to locate the page associated with the amino acid name on the PDBE site.

    Returns
    -------
    aa_dict : dict
        Keys are AminoAcidDB field names.
        Values are the str values for the new amino acid, scraped from the PDBE if possible. None if not found.
    """
    letter = 'X'
    # Try to get content from PDBE.
    url_string = 'http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{0}'.format(code)
    r = requests.get(url_string)
    # Raise error if content not obtained.
    if not r.ok:
        raise IOError("Could not get to url {0}".format(url_string))

    # Parse r.text in an ugly way to get the required information.
    description = r.text.split('<h3>Molecule name')[1].split('</tr>')[0]
    description = description.strip().split('\n')[3].strip()[:255]
    modified = r.text.split("<h3>Standard parent ")[1].split('</tr>')[0]
    modified = modified.replace(" ", "").replace('\n', '').split('<')[-3].split('>')[-1]
    if modified == "NotAssigned":
        modified = None
    # Add the required information to a dictionary which can then be passed to add_amino_acid_to_json.
    aa_dict = {'code': code, 'description': description, 'modified': modified, 'letter': letter}
    return aa_dict


def add_amino_acid_to_json(code, description, letter='X', modified=None, force_add=False):
    """ Add an amino acid to the amino_acids.json file used to populate the amino_acid table.

    Parameters
    ----------
    code : str
        New code to be added to amino acid table.
    description : str
        Description of the amino acid, e.g. 'amidated terminal carboxy group'.
    letter : str, optional
        One letter code for the amino acid.
        Defaults to 'X'
    modified : str or None, optional
        Code of modified amino acid, e.g. 'ALA', or None.
        Defaults to None
    force_add : bool, optional
        If True, will over-write existing dictionary value for code if already in amino_acids.json.
        If False, then an IOError is raised if code is already in amino_acids.json.

    Raises
    ------
    IOError
        If code is already in amino_acids.json and force_add is False.

    Returns
    -------
    None

    """
    # If code is already in the dictionary, raise an error
    if (not force_add) and code in amino_acids_dict.keys():
        raise IOError("{0} is already in the amino_acids dictionary, with values: {1}".format(
            code, amino_acids_dict[code]))

    # Prepare data to be added.
    add_code = code
    add_code_dict = {'description': description, 'letter': letter, 'modified': modified}

    # Check that data does not already exist, and if not, add it to the dictionary.
    amino_acids_dict[add_code] = add_code_dict

    # Write over json file with updated dictionary.
    with open(_amino_acids_json_path, 'w') as foo:
        foo.write(json.dumps(amino_acids_dict))

    return


__author__ = 'Jack W. Heal'
