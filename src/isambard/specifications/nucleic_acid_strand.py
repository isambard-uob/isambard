"""Contains code for generating single polynucleotides."""

from collections import OrderedDict

import numpy
from ampal.geometry import (angle_between_vectors, unit_vector, find_foot,
                            Quaternion, dihedral, find_transformations,
                            cylindrical_to_cartesian, Axis)
from ampal import Atom, Polynucleotide, Nucleotide

_helix_parameters = {
    # nucleotides_per_turn, rise_per_nucleotide
    'a_dna': (11.0, 2.6),
    'b_dna': (10.2, 3.4),
    'z_dna': (12.3, 3.7)
}


_helix_handedness = {
    'a_dna': 'r',
    'b_dna': 'r',
    'z_dna': 'r'
}


_backbone_properties = {
    # helix_type : {atom_type: [radius, angular_offset (in radians), z-shift]}
    'a_dna': {
        'P': [8.92, 0, 0],  # P
    },
    'b_dna': {
        'mol_code_format': 'D{}',
        'labels': ["P", "OP1", "OP2", "O5'", "C5'", "C4'",
                   "O4'", "C3'", "O3'", "C2'", "C1'"],
        'atoms': [
            (8.92, 0.0, 0.0),
            (10.23, 0.07, 0.21),
            (9.0, -0.14, 0.77),
            (7.72, -0.5, -3.11),
            (7.68, -0.36, -4.05),
            (7.57, -0.17, -3.33),
            (6.22, -0.11, -3.17),
            (8.19, -0.19, -1.93),
            (8.74, -0.04, -1.55),
            (7.02, -0.25, -1.06),
            (5.79, -0.14, -1.83),
            (4.58, -0.30, -1.83)
        ]
    },
    'z_dna': {
        'P': [8.92, 0, 0],  # P
    }
}

_bases = {
    'G': {
        'labels': ["O4'", "C1'", "N9", "C8", "N7", "C5",
                   "C6", "O6", "N1", "C2", "N2", "N3", "C4"],
        'ref_atom': "N9",
        'rot_adj': 15.0,
        'atoms': [
            (33.930, 30.325, 17.814),
            (34.631, 29.403, 17.056),
            (35.908, 29.190, 17.720),
            (36.117, 28.636, 18.964),
            (37.376, 28.606, 19.309),
            (38.040, 29.196, 18.238),
            (39.419, 29.435, 18.041),
            (40.348, 29.171, 18.808),
            (39.672, 30.068, 16.810),
            (38.706, 30.393, 15.889),
            (39.131, 30.975, 14.756),
            (37.406, 30.170, 16.069),
            (37.150, 29.567, 17.255)
        ]
    },
    'C': {
        'labels': ["O4'", "C1'", "N1", "C2", "O2", "N3",
                   "C4", "N4", "C5", "C6"],
        'ref_atom': "N1",
        'rot_adj': 47.0,
        'atoms': [
            (35.422, 27.992, 13.920),
            (36.655, 27.383, 13.630),
            (37.240, 26.863, 14.903),
            (38.622, 26.960, 15.116),
            (39.337, 27.469, 14.237),
            (39.137, 26.506, 16.285),
            (38.336, 25.972, 17.205),
            (38.892, 25.545, 18.338),
            (36.925, 25.864, 17.008),
            (36.426, 26.317, 15.853)
        ]
    },
    'A': {
        'labels': ["O4'", "C1'", "N9", "C8", "N7", "C5",
                   "C6", "N6", "N1", "C2", "N3", "C4"],
        'ref_atom': "N9",
        'rot_adj': 5.0,
        'atoms': [
            (39.099, 26.091,  9.763),
            (39.785, 24.859,  9.897),
            (39.558, 24.371, 11.250),
            (38.407, 23.814, 11.737),
            (38.478, 23.462, 12.998),
            (39.769, 23.803, 13.357),
            (40.461, 23.687, 14.561),
            (39.913, 23.168, 15.662),
            (41.745, 24.130, 14.598),
            (42.280, 24.652, 13.481),
            (41.715, 24.823, 12.282),
            (40.447, 24.374, 12.291)
        ]
    },
    'T': {
        'labels': ["O4'", "C1'", "N1", "C2", "O2", "N3",
                   "C4", "O4", "C5", "C7", "C6"],
        'ref_atom': "N1",
        'rot_adj': 30.0,
        'atoms': [
            (43.43,  23.953,  9.349),
            (44.116, 22.801,  9.767),
            (43.234, 22.076, 10.744),
            (43.702, 21.817, 12.016),
            (44.824, 22.102, 12.389),
            (42.814, 21.185, 12.840),
            (41.527, 20.804, 12.542),
            (40.808, 20.241, 13.364),
            (41.084, 21.115, 11.196),
            (39.697, 20.746, 10.751),
            (41.948, 21.737, 10.374)
        ]
    }
}


class NucleicAcidStrand(Polynucleotide):
    """Generates a `Polynucleotide` from a sequence of bases.

    Parameters
    ----------
    sequence: str
        The nucleotide sequence of the nucleic acid.
    helix_type: str, optional
        The type of nucleic acid helix to generate.
    phos_3_prime: bool, optional
        If false the 5' and the 3' phosphor will be omitted.

    Attributes
    ----------
    base_sequence : str
        Nucleotide sequence for the nucleic acid.
    num_monomers : int
        Number of bases in the sequence.
    helix_type : str
        The type of nucleic acid helix. Currently only b-DNA is
        implemented.
    phos_3_prime : bool
        If `true`, the 3' end of the strand will have a phosphate.
    nucleotides_per_turn : float
        Number of nucleotides per turn of the nucleic acid super
        helix.
    rise_per_nucleotide : float
        Rise along the nucleic acid super helix per base.
    handedness : str
        Handedness of the super helix.
    helix_start : 3D Vector (tuple or list or numpy.array)
        Initial coordinate of the backbone primitive.
    helix_end : 3D Vector (tuple or list or numpy.array)
        Last coordinate of the backbone primitive.

    Raises
    ------
    ValueError
        Raised if a sequence is not provided.
    """

    def __init__(self, sequence='GATC', helix_type='b_dna', phos_3_prime=False):
        super().__init__()
        if not sequence:
            raise ValueError(
                'A sequence must be provided to build a region of DNA.')
        self.base_sequence = sequence
        self.num_monomers = len(sequence)
        self.helix_type = helix_type
        self.phos_3_prime = phos_3_prime
        self.nucleotides_per_turn, self.rise_per_nucleotide = _helix_parameters[
            self.helix_type]
        self.handedness = _helix_handedness[self.helix_type]
        self.helix_start = numpy.array([0.0, 0.0, 0.0])
        self.helix_end = self.helix_length * numpy.array([0.0, 0.0, 1.0])
        self.build()

    @classmethod
    def from_start_and_end(cls, start, end, sequence, helix_type='b_dna',
                           phos_3_prime=False):
        """Generates a helical `Polynucleotide` that is built along an axis.

        Parameters
        ----------
        start: [float, float, float]
            Start of the build axis.
        end: [float, float, float]
            End of build axis.
        sequence: str
            The nucleotide sequence of the nucleic acid.
        helix_type: str
            The type of nucleic acid helix to generate.
        phos_3_prime: bool
            If false the 5' and the 3' phosphor will be omitted.
        """
        start = numpy.array(start)
        end = numpy.array(end)
        instance = cls(sequence, helix_type=helix_type,
                       phos_3_prime=phos_3_prime)
        instance.move_to(start=start, end=end)
        return instance

    @property
    def axis(self):
        """The super-helical axis."""
        return Axis(start=self.helix_start, end=self.helix_end)

    @property
    def ax_unit(self):
        """The unit tangent of the super-helical axis."""
        return self.axis.unit_tangent

    @property
    def rad_unit(self):
        """The unit normal of the super-helical axis."""
        return self.axis.unit_normal

    @property
    def tan_unit(self):
        """The unit binormal of the super-helical axis."""
        return self.axis.unit_binormal

    @property
    def helix_length(self):
        """Length of the helix in Ångstroms."""
        return self.num_monomers * self.rise_per_nucleotide

    def translate(self, vector):
        super().translate(vector=vector)
        self.helix_start += vector
        self.helix_end += vector
        return

    def rotate(self, angle, axis, point=None, radians=False):
        super().rotate(
            angle=angle, axis=axis, point=point, radians=radians)
        # modify helix_start and helix_end accordingly.
        q = Quaternion.angle_and_axis(angle=angle, axis=axis, radians=radians)
        self.helix_start = q.rotate_vector(v=self.helix_start, point=point)
        self.helix_end = q.rotate_vector(v=self.helix_end, point=point)
        return

    def build(self):
        """Build single DNA strand along z-axis, starting with P on x-axis"""
        ang_per_res = (2 * numpy.pi) / self.nucleotides_per_turn
        atom_offset_coords = _backbone_properties[self.helix_type]['atoms']
        if self.handedness == 'l':
            handedness = -1
        else:
            handedness = 1
        base_atom_labels = _backbone_properties[self.helix_type]['labels']
        monomers = []
        mol_code_format = _backbone_properties[self.helix_type]['mol_code_format']
        for i, b in enumerate(self.base_sequence):
            nucleotide = Nucleotide(
                mol_code=mol_code_format.format(b), parent=self)
            atoms_dict = OrderedDict()
            if (i == (len(self.base_sequence) - 1)) and not self.phos_3_prime:
                # Do not include phosphate on last nucleotide
                atom_labels = base_atom_labels[3:] + [_bases[b]['labels'][2]]
                atom_offsets = {k: v for k, v in zip(
                    atom_labels, atom_offset_coords[3:])}
            else:
                atom_labels = base_atom_labels + [_bases[b]['labels'][2]]
                atom_offsets = {k: v for k, v in zip(
                    atom_labels, atom_offset_coords)}
            for atom_label in atom_labels:
                r, zeta, z_shift = atom_offsets[atom_label]
                rot_ang = ((i * ang_per_res) + zeta) * handedness
                z = (self.rise_per_nucleotide * i) + z_shift
                coords = cylindrical_to_cartesian(
                    radius=r, azimuth=rot_ang, z=z, radians=True)
                atom = Atom(
                    coordinates=coords, element=atom_label[0],
                    parent=nucleotide, res_label=atom_label)
                atoms_dict[atom_label] = atom
            base_ref = _bases[b]['ref_atom']
            rot_adj = _bases[b]['rot_adj']
            base_dict = OrderedDict(
                zip(_bases[b]['labels'], _bases[b]['atoms']))
            translation, angle, axis, point = find_transformations(
                base_dict[base_ref], base_dict["C1'"],
                atoms_dict[base_ref]._vector, atoms_dict["C1'"]._vector)
            q1 = Quaternion.angle_and_axis(angle, axis)
            # Align N9 C1'
            for k, v in base_dict.items():
                base_dict[k] = q1.rotate_vector(v, point) + translation
            # Rotate to align O4'
            axis = numpy.array(base_dict["C1'"]) - base_dict[base_ref]
            angle = dihedral(base_dict["O4'"], base_dict[base_ref],
                             base_dict["C1'"], atoms_dict["O4'"]) - rot_adj
            q2 = Quaternion.angle_and_axis(angle, axis)
            for k, v in list(base_dict.items()):
                if k not in atoms_dict:
                    atom = Atom(q2.rotate_vector(v, base_dict[base_ref]),
                                element=k[0], parent=nucleotide,
                                res_label=k)
                    atoms_dict[k] = atom
            nucleotide.atoms = atoms_dict
            monomers.append(nucleotide)
        self._monomers = monomers
        self.relabel_monomers()
        self.relabel_atoms()
        return

    def move_to(self, start, end):
        """Moves the `Polynucleotide` to lie on the `start` and `end` vector.

        Parameters
        ----------
        start : 3D Vector (tuple or list or numpy.array)
            The coordinate of the start of the helix primitive.
        end : 3D Vector (tuple or list or numpy.array)
            The coordinate of the end of the helix primitive.

        Raises
        ------
        ValueError
            Raised if `start` and `end` are very close together.
        """
        start = numpy.array(start)
        end = numpy.array(end)
        if numpy.allclose(start, end):
            raise ValueError('start and end must NOT be identical')
        translation, angle, axis, point = find_transformations(
            self.helix_start, self.helix_end, start, end)
        if not numpy.isclose(angle, 0.0):
            self.rotate(angle=angle, axis=axis, point=point, radians=False)
        self.translate(vector=translation)
        return


__author__ = 'Jack W. Heal, Christopher W. Wood'
