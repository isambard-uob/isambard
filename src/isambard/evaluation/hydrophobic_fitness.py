"""Module for calculating the hydrophobic fitness of a protein."""


import ampal
from ampal.amino_acids import standard_amino_acids


HYDROPHOBIC = ['C', 'F', 'I', 'L', 'M', 'V', 'W']


def calculate_hydrophobic_fitness(assembly):
    r"""Calculates the hydrophobic fitness of a protein.

    Hydrophobic fitness is an efficient centroid-based method for
    calculating the packing quality of your structure [3]_. For this
    method C, F, I, L, M, V, W and Y are considered hydrophobic. The
    algorithm has two terms:

    .. math::

        Hydrophobic\ term = \frac{\sum\limits_{i}^{}(H_{i}-H_{i}^{\circ})}{n}

    where :math:`H_{i}` is the total number of number of hydrophobic
    contacts of residue *i*, :math:`H_{i}^{\circ}` is the number of
    hydrophobic contacts expected by chance :eq:`hi_chance` and *n*
    is the total number of residues.

    .. math::

        Burial\ term = \frac{\sum\limits_{i}B_{i}}{n}

    where :math:`B_{i}` is the number of centroids within 10 A. The
    number of hydrophobic contacts expected by chance is calculated
    as follows:

    .. math::
        :label: hi_chance

        H_{i}^{\circ} = C_{i}\left(\frac{h_{i}}{N_{i}}\right)

    where :math:`C_{i}` is the number of all side-chain centroids
    that contact residue *i*, *h* is the total number of hydrophobic
    residues in the sequence except for any neighbours. :math:`h_{i}`
    is the total number of residue minus i and neighbours. The
    hydrophobic fitness score is the combination of these terms:

    .. math::

        HF = -\frac{\left(\sum\limits_{i}B_{i}\right)\left(\sum\limits_{
        i}^{}(H_{i}-H_{i}^{\circ})\right)}{n^{2}}

    Notes
    -----
    WARNING: The scores produced by this implementation do not quite
    match scores listed in publications from the Levitt group. The
    scores are generally off by up to around 1 unit, and so it
    should still be useful.

    References
    ----------
    .. [3] Huang ES, Subbiah S and Levitt M (1995) Recognizing Native
       Folds by the Arrangement of Hydrophobic and Polar Residues, J. Mol.
       Biol return., **252**, 709-720.
    """
    hydrophobic_centroids = []
    tyrosine_centroids = []
    polar_centroids = []
    for residue in [r for r in assembly.get_monomers()
                    if isinstance(r, ampal.Residue)]:
        centroid_list = None
        centroid = residue.centroid
        if residue.mol_letter in HYDROPHOBIC:
            centroid_list = hydrophobic_centroids
        elif residue.mol_letter == 'Y':
            centroid_list = tyrosine_centroids
        elif residue.mol_letter in standard_amino_acids:
            centroid_list = polar_centroids
        else:
            continue
        if centroid_list is not None:
            centroid_list.append(
                (residue.parent.id, int(residue.id),
                 residue['CA'] if centroid is None else centroid))
    hf = run_hf_loop(hydrophobic_centroids,
                     tyrosine_centroids, polar_centroids)
    return hf


def run_hf_loop(hydrophobic_centroids, tyrosine_centroids, polar_centroids):
    """Runs the hydrophobic fitness algorithm.

    Parameters
    ----------
    hydrophobic_centroids : [(str, int, (float, float, float))]
        A list containing the chain ID, residue number and the centroid
        coordinate position for all the hydrophobic residues exclusing
        tyrosine.
    tyrosine_centroids : [(str, int, (float, float, float))]
        A list containing the chain ID, residue number and the centroid
        coordinate position for all the tyrosine residues.
    polar_centroids : [(str, int, (float, float, float))]
        A list containing the chain ID, residue number and the centroid
        coordinate position for all the polar residues.

    Returns
    -------
    hydrophobic_fitness : float
        The hydrophobic fitness score.
    """
    n_hydrophobic = len(hydrophobic_centroids)
    n_tyrosine = len(tyrosine_centroids)
    n_polar = len(polar_centroids)
    hydrophobic_scores = []
    burial_scores = []
    for reference in hydrophobic_centroids:
        h_centroids_in_7_3, h_centroids_in_10, h_neighbours = get_number_within(
            reference, hydrophobic_centroids)
        y_centroids_in_7_3, y_centroids_in_10, y_neighbours = get_number_within(
            reference, tyrosine_centroids)
        p_centroids_in_7_3, p_centroids_in_10, p_neighbours = get_number_within(
            reference, polar_centroids)
        all_neighbours = h_neighbours + y_neighbours + p_neighbours
        ci = (h_centroids_in_7_3 + y_centroids_in_7_3 + p_centroids_in_7_3) - (
            all_neighbours)
        hi = (n_hydrophobic - h_neighbours) + (n_tyrosine - y_neighbours)
        Ni = hi + (n_polar - p_neighbours)
        Hic = ci * (hi / Ni)
        Hi = (h_centroids_in_7_3 + y_centroids_in_7_3) - (
            h_neighbours + y_neighbours)
        hydrophobic_scores.append(Hi - Hic)
        burial_scores.append(
            (h_centroids_in_10 + y_centroids_in_10 + p_centroids_in_10))
    hydrophobic_fitness = -1 * (
        (sum(burial_scores)*sum(hydrophobic_scores)
         ) / ((n_hydrophobic+n_tyrosine)**2)
    )
    return hydrophobic_fitness


def get_number_within(reference, target_points):
    """Get the number of points within 10 and 7.3 residue.

    Parameters
    ----------
    reference : (str, int, (float, float, float))
        Reference centroid.
    target_points : [(str, int, (float, float, float))]
        A list of centroids.

    Returns
    -------
    within_and_neighbours : (int, int, int)
        Returns the number of target centroids within 7.3 A, 10.0 A
        and the number of neighbours (based on chain ID and residue
        number).
    """
    centroids_in_10 = 0
    centroids_in_7_3 = 0
    neighbours = 0
    for target in target_points:
        if target is reference:
            continue
        elif reference[0] == target[0]:
            if ((reference[1] - 1) == target[1]) or (
                    (reference[1] + 1) == target[1]):
                neighbours += 1
        distance = ampal.geometry.distance(reference[2], target[2])
        if distance <= 10.0:
            centroids_in_10 += 1
        if distance <= 7.3:
            centroids_in_7_3 += 1
    return centroids_in_7_3, centroids_in_10, neighbours


__author__ = 'Christopher W. Wood'
