"""Contains functionality for builting protein structure from torsion angles."""
#cython: embedsignature=True

from collections import OrderedDict
import numpy
from libc.math cimport cos, sin, M_PI

from ampal import Atom, Polypeptide, Residue, align
from ampal.analyse_protein import measure_torsion_angles
from ampal.amino_acids import ideal_backbone_bond_lengths, ideal_backbone_bond_angles
from ampal.geometry import dihedral, Quaternion


cdef _planar_c_and_o_coords(double n_ca, ca_c, c_o, n_ca_c, ca_c_o):
    cdef double cy, cz, oy, oz, half_pi
    n_ca_c = (n_ca_c / 180.0) * M_PI
    ca_c_o = (ca_c_o / 180.0) * M_PI
    half_pi = M_PI / 2
    cy = cos(n_ca_c - half_pi) * ca_c
    cz = sin(n_ca_c - half_pi) * ca_c + n_ca
    oy = cy + cos(half_pi + n_ca_c - ca_c_o) * c_o
    oz = cz + sin(half_pi + n_ca_c - ca_c_o) * c_o
    return cy, cz, oy, oz


def planar_residue(n_ca=1.47, ca_c=1.53, c_o=1.24, n_ca_c=110.0, ca_c_o=121.0):
    """
    planar_residue(n_ca=1.47, ca_c=1.53, c_o=1.24, n_ca_c=110.0, ca_c_o=121.0)

    GLY Residue with all coordinates in the y-z plane.

    Note
    ----
    Builds a GLY Residue in the y-z plane.
    N is placed at the origin, and CA on the positive z-axis.
    The C and O coordinates are determined by input bond length and angles.

    Parameters
    ----------
    n_ca: float, optional
        N-CA bond distance in Angstroms.
    ca_c: float, optional
        CA-C bond distance in Angstroms.
    c_o: float, optional
        C=O bond distances in Angstroms.
    n_ca_c: float, optional
        N-CA-C bond angle in degrees.
    ca_c_o: float, optional
        CA-C=O bond angle in degrees.

    Returns
    -------
    Residue
        A GLY Residue in the y-z plane.
    """
    n = Atom([0.0, 0.0, 0.0], element='N', res_label='N')
    ca = Atom([0.0, 0.0, n_ca], element='C', res_label='CA')
    cy, cz, oy, oz = _planar_c_and_o_coords(
        n_ca=n_ca, ca_c=ca_c, c_o=c_o, n_ca_c=n_ca_c, ca_c_o=ca_c_o)
    c = Atom([0.0, cy, cz], element='C', res_label='C')
    o = Atom(numpy.array([0.0, oy, oz]), element='O', res_label='O')
    atoms = OrderedDict([('N', n), ('CA', ca), ('C', c), ('O', o)])
    res = Residue(atoms=atoms, mol_code='GLY')
    for atom in res.get_atoms():
        atom.ampal_parent = res
    return res


class TAPolypeptide(Polypeptide):
    """Generates a polypeptide with defined torsion angles.

    Note
    ----
    Default bond angles and distances from [1]. The final oxygen atom
    is positioned as if the next (imaginary) Residue is bonded in a
    trans conformation (i.e. using a default value for the CA-C=O 
    angle of 121.0, not 119.0).

    References
    ----------
    .. [0] Schulz, G. E, and R. Heiner Schirmer. Principles Of Protein
       Structure. New York: Springer-Verlag, 1979.

    Parameters
    ----------
    torsion_angles : list of triples
        Contains lists of torsion angles (omega, phi, psi) that will
        be used to create the polypeptide. The number of residues
        is equal to the length of this list. The first value of the
        first list (omega) and the last value of the last list (psi)
        will not be used to generate the polypeptide.
    auto_build : bool
        If `true`, the model will be built as part of instantiation.

    Attributes
    ----------
    torsion_angles : list of triples
        Contains lists of torsion angles (omega, phi, psi) that will
        be used to create the polypeptide. The number of residues
        is equal to the length of this list. The first value of the
        first list (omega) and the last value of the last list (psi)
        will not be used to generate the polypeptide.
    n_ca_bonds: list of float
        Backbone N-CA bond distances in Angstroms.
    ca_c_bonds: list of float
        Backbone CA-C bond distances in Angstroms.
    c_n_bonds: list of float
        Backbone C-N (peptide) bond distances in Angstroms.
    c_o_bonds: list of float
        Backbone C=O bond distances in Angstroms.
    n_ca_c_angles: list of float
        Backbone N-CA-C bond angle in degrees.
    ca_c_o_angles: list of float
        Backbone CA-C=O bond angle in degrees.
    ca_c_n_angles: list of float
        Backbone CA-C-N bond angle in degrees.
    c_n_ca_angles: list of float
        Backbone C-N-CA bond angle in degrees.
    """

    def __init__(self, torsion_angles, auto_build=True):
        super().__init__()
        self.torsion_angles = torsion_angles
        self.num_monomers = len(self.torsion_angles)
        self.trans = [True if abs(self.torsion_angles[i][0]) >= 90
                      else False
                      for i in range(self.num_monomers)]
        # backbone bond lengths
        self.c_o_bonds = [ideal_backbone_bond_lengths['c_o']
                          for _ in range(self.num_monomers)]
        self.ca_c_bonds = [ideal_backbone_bond_lengths['ca_c']
                           for _ in range(self.num_monomers)]
        self.n_ca_bonds = [ideal_backbone_bond_lengths['n_ca']
                           for _ in range(self.num_monomers)]
        self.c_n_bonds = [ideal_backbone_bond_lengths['c_n']
                          for _ in range(self.num_monomers - 1)]
        # backbone bond angles
        self.n_ca_c_angles = [
            ideal_backbone_bond_angles['trans']['n_ca_c'] if self.trans[i]
            else ideal_backbone_bond_angles['cis']['n_ca_c']
            for i in range(self.num_monomers)]
        self.ca_c_o_angles = [
            ideal_backbone_bond_angles['trans']['ca_c_o'] if self.trans[i + 1]
            else ideal_backbone_bond_angles['cis']['ca_c_o']
            for i in range(self.num_monomers - 1)]
        self.ca_c_o_angles.append(
            ideal_backbone_bond_angles['trans']['ca_c_o'])
        self.ca_c_n_angles = [
            ideal_backbone_bond_angles['trans']['ca_c_n'] if self.trans[i + 1]
            else ideal_backbone_bond_angles['cis']['ca_c_n']
            for i in range(self.num_monomers - 1)]
        self.c_n_ca_angles = [
            ideal_backbone_bond_angles['trans']['c_n_ca'] if self.trans[i + 1]
            else ideal_backbone_bond_angles['cis']['c_n_ca']
            for i in range(self.num_monomers - 1)]
        if auto_build:
            self.build()

    @classmethod
    def from_polypeptide(cls, polypeptide, align_model=True, auto_build=True):
        """ Builds input polypeptide as a TAPolypeptide

        Parameters
        ----------
        polypeptide : Polypeptide
        align_model : bool
            If True, runs align() on the TAPolypeptide to align it with
            the Polypeptide before returning.
        auto_build : bool
            If `true`, the model will be built as part of instantiation.

        Raises
        ------
        TypeError
            Raised if `polypeptide` is not `ampal.Polypeptide`.
        """
        if not isinstance(polypeptide, Polypeptide):
            raise TypeError('Invalid object type {0} given. This method can '
                            'only be used with a Polypeptide'
                            .format(type(polypeptide)))
        torsion_angles = [list(x) for x in measure_torsion_angles(polypeptide)]
        # dummy value for first omega (not used in build, but can't be None)
        torsion_angles[0][0] = 180.0
        instance = cls(torsion_angles=torsion_angles, auto_build=False)
        bond_angles = polypeptide.backbone_bond_angles
        bond_lengths = polypeptide.backbone_bond_lengths
        instance.n_ca_bonds = bond_lengths['n_ca']
        instance.ca_c_bonds = bond_lengths['ca_c']
        instance.c_o_bonds = bond_lengths['c_o']
        instance.c_n_bonds = bond_lengths['c_n']
        instance.n_ca_c_angles = bond_angles['n_ca_c']
        instance.ca_c_o_angles = bond_angles['ca_c_o']
        instance.ca_c_n_angles = bond_angles['ca_c_n']
        instance.c_n_ca_angles = bond_angles['c_n_ca']
        if auto_build:
            instance.build()
            # Rotate about CA-C bond to align O atoms
            # (there is some rotational freedom of this atom in natural structures).
            polypeptide_n_ca_c_o_angles = [
                dihedral(r['N'], r['CA'], r['C'], r['O'])
                for r in polypeptide.get_monomers(ligands=False)]
            instance_n_ca_c_o_angles = [
                dihedral(r['N'], r['CA'], r['C'], r['O'])
                for r in instance.get_monomers()]
            for i, r in enumerate(instance):
                angle = (polypeptide_n_ca_c_o_angles[i]
                         - instance_n_ca_c_o_angles[i])
                q = Quaternion.angle_and_axis(
                    angle=angle, axis=(r['C'] - r['CA']))
                r['O']._vector = q.rotate_vector(
                    r['O']._vector, point=r['C']._vector)
            if align_model:
                align(target=polypeptide, mobile=instance)
        return instance

    @property
    def o_c_n_angles(self):
        """Backbone O=C-N angle determined by CA-C=O and CA-C-N angle.

        Notes
        -----
        Must sum to 360.0.
        """
        return [360.0 - (x + y)
                for x, y in zip(self.ca_c_o_angles, self.ca_c_n_angles)]

    def build(self):
        """Builds TAPolypeptide using its bond angle and bond length lists."""
        planar_residues = [planar_residue(
            n_ca=self.n_ca_bonds[i],
            ca_c=self.ca_c_bonds[i],
            c_o=self.c_o_bonds[i],
            n_ca_c=self.n_ca_c_angles[i],
            ca_c_o=self.ca_c_o_angles[i])
            for i in range(self.num_monomers)]
        r = Polypeptide([planar_residues[0]])
        for i in range(self.num_monomers - 1):
            r2 = planar_residues[i + 1]
            omega, phi = self.torsion_angles[i + 1][:2]
            psi = self.torsion_angles[i][2]
            r.c_join(r2, omega=omega, phi=phi, psi=psi,
                     o_c_n_angle=self.o_c_n_angles[i],
                     c_n_ca_angle=self.c_n_ca_angles[i],
                     c_n_length=self.c_n_bonds[i], relabel=False)
        r.relabel_all()
        self._monomers = r._monomers
        for monomer in self._monomers:
            monomer.ampal_parent = self
        return


__author__ = 'Jack W. Heal, Christopher W. Wood'
