"""Contains classes and functions for modelling helical structure."""

from collections import OrderedDict

import numpy
from ampal.geometry import (Quaternion, dihedral, find_transformations,
                            cylindrical_to_cartesian, Axis, HelicalCurve)
from ampal import Atom, Residue, Polypeptide
from ampal.pseudo_atoms import Primitive

_helix_parameters = {
    # residues_per_turn, rise_per_residue
    'alpha': (3.613, 1.52),
    'pi': (4.385, 1.202),
    '3-10': (3.0, 2.0),
    'PPI': (3.3, 1.9),
    'PPII': (3.0, 3.1),
    'collagen': (3.3, 2.9),
}


_helix_handedness = {
    'alpha': 'r',
    'pi': 'r',
    '3-10': 'r',
    'PPI': 'r',
    'PPII': 'l',
    'collagen': 'l',
}


_atom_offsets = {
    # helix_type : {atom_type: [radius, angular_offset (in radians), z-shift]}
    'alpha': {
        'N': [1.50, -0.47, -0.90],  # N
        'CA': [2.30, 0.00, 0.00],  # CA
        'C': [1.71, 0.48, 1.05],  # C
        'O': [2.08, 0.42, 2.22],  # O
        #'CB': [3.32, 0.27, -0.86],  # CB
    },
    'pi': {
        'N': [2.014, -0.433, -0.725],  # N
        'CA': [2.786, 0.000, 0.000],  # CA
        'C': [2.242, 0.424, 0.966],  # C
        'O': [2.530, 0.371, 2.166],  # O
        #'CB': [3.715, 0.219, -0.977],  # CB
    },
    '3-10': {
        'CA': [2.30, 0.00, 0.00],  # CA
    },
    'PPI': {
        'CA': [2.30, 0.00, 0.00],  # CA
    },
    'PPII': {
        'N': [1.043, -0.732, -1.172],  # N
        'CA': [1.327, 0.000, 0.000],  # CA
        'C': [0.358, -0.265, 1.170],  # C
        'O': [1.178, -1.793, 1.413],  # O
        #'CB': [2.741, -0.214, -0.383], # CB
    },
    'collagen': {
        'N': [1.043, -0.732, -1.172],  # N
        'CA': [1.327, 0.000, 0.000],  # CA
        'C': [0.358, -0.265, 1.170],  # C
        'O': [1.178, -1.793, 1.413],  # O
        # 'CB': [2.741, -0.214, -0.383], # CB
    },
}


class Helix(Polypeptide):
    """Creates a model of a `Polypeptide` with helical structure.

    Parameters
    ----------
    aa : int, optional
        Number of amino acids in the `Helix`.
    helix_type : str, optional
        Type of helix, can be: 'alpha', 'pi', '3-10',
        'PPI', 'PPII', 'collagen'.

    Attributes
    ----------
    num_monomers : int
        Number of amino acids in the `Helix`.
    helix_type : str
        Type of helix, can be: 'alpha', 'pi', '3-10',
        'PPI', 'PPII', 'collagen'.
    residues_per_turn : float
        Number of residues per turn of the helix.
    rise_per_residue : float
        Distance that each residue travels along the z-axis.
    handedness : str
        The handedness of the helix can be 'r' or 'l'.
    helix_start : 3D Vector (tuple or list or numpy.array)
        The coordinate of the start of the helix primitive.
    helix_end : 3D Vector (tuple or list or numpy.array)
        The coordinate of the end of the helix primitive.
    """

    def __init__(self, aa=10, helix_type='alpha'):
        super(Helix, self).__init__()
        self.num_monomers = aa
        self.helix_type = helix_type
        self.residues_per_turn, self.rise_per_residue = _helix_parameters[self.helix_type]
        self.handedness = _helix_handedness[self.helix_type]
        self.helix_start = numpy.array([0.0, 0.0, 0.0])
        self.helix_end = self.helix_length * numpy.array([0.0, 0.0, 1.0])
        self.build()

    def __repr__(self):
        if len(self.sequence) > 15:
            seq = self.sequence[:12] + '...'
        else:
            seq = self.sequence
        return '<Helix containing {} {}. Sequence: {}>'.format(
            len(self._monomers), 'Residue' if len(self._monomers) == 1 else 'Residues', seq)

    @classmethod
    def from_start_and_end(cls, start, end, aa=None, helix_type='alpha'):
        """Creates a `Helix` between `start` and `end`.

        Parameters
        ----------
        start : 3D Vector (tuple or list or numpy.array)
            The coordinate of the start of the helix primitive.
        end : 3D Vector (tuple or list or numpy.array)
            The coordinate of the end of the helix primitive.
        aa : int, optional
            Number of amino acids in the `Helix`. If `None, an
            appropriate number of residues are added.
        helix_type : str, optional
            Type of helix, can be: 'alpha', 'pi', '3-10',
            'PPI', 'PPII', 'collagen'.
        """
        start = numpy.array(start)
        end = numpy.array(end)
        if aa is None:
            rise_per_residue = _helix_parameters[helix_type][1]
            aa = int((numpy.linalg.norm(end - start) / rise_per_residue) + 1)
        instance = cls(aa=aa, helix_type=helix_type)
        instance.move_to(start=start, end=end)
        return instance

    @property
    def axis(self):
        """The `Axis` around which the helix is constructed."""
        return Axis(start=self.helix_start, end=self.helix_end)

    @property
    def ax_unit(self):
        """The unit tangent of the axis."""
        return self.axis.unit_tangent

    @property
    def rad_unit(self):
        """The unit normal of the axis."""
        return self.axis.unit_normal

    @property
    def tan_unit(self):
        """The unit binormal of the axis."""
        return self.axis.unit_binormal

    @property
    def helix_length(self):
        """The length of the `Helix`."""
        return self.num_monomers * self.rise_per_residue

    def translate(self, vector):
        super(Helix, self).translate(vector=vector)
        self.helix_start += vector
        self.helix_end += vector
        return

    def rotate(self, angle, axis, point=None, radians=False):
        super(Helix, self).rotate(angle=angle,
                                  axis=axis, point=point, radians=radians)
        # modify helix_start and helix_end accordingly.
        q = Quaternion.angle_and_axis(angle=angle, axis=axis, radians=radians)
        self.helix_start = q.rotate_vector(v=self.helix_start, point=point)
        self.helix_end = q.rotate_vector(v=self.helix_end, point=point)
        return

    # TODO Separate out the polypeptide bit from the geometry.
    def build(self):
        """Build straight helix along z-axis, starting with CA1 on x-axis"""
        ang_per_res = (2 * numpy.pi) / self.residues_per_turn
        atom_offsets = _atom_offsets[self.helix_type]
        if self.handedness == 'l':
            handedness = -1
        else:
            handedness = 1

        atom_labels = ['N', 'CA', 'C', 'O']
        if all([x in atom_offsets.keys() for x in atom_labels]):
            res_label = 'GLY'
        else:
            res_label = 'UNK'

        monomers = []
        for i in range(self.num_monomers):
            residue = Residue(mol_code=res_label, parent=self)
            atoms_dict = OrderedDict()
            for atom_label in atom_labels:
                r, zeta, z_shift = atom_offsets[atom_label]
                rot_ang = ((i * ang_per_res) + zeta) * handedness
                z = (self.rise_per_residue * i) + z_shift
                coords = cylindrical_to_cartesian(
                    radius=r, azimuth=rot_ang, z=z, radians=True)
                atom = Atom(
                    coordinates=coords, element=atom_label[0],
                    parent=residue, res_label=atom_label)
                atoms_dict[atom_label] = atom
            residue.atoms = atoms_dict
            monomers.append(residue)
        self._monomers = monomers
        self.relabel_monomers()
        self.relabel_atoms()
        return

    def move_to(self, start, end):
        """Moves the `Helix` to lie on the vector between `start` and `end`.

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


class HelicalHelix(Polypeptide):
    """Builds a `Helix` along a curve defined by a super helix.

    Parameters
    ----------
    aa : int
        Number of amino acids in the `Helix`. If `None, an
        appropriate number of residues are added.
    major_pitch : float, optional
        Pitch of the super helix.
    major_radius : float, optional
        Radius of the super helix.
    major_handedness : str, optional
        Handedness of the super helix.
    minor_helix_type : str, optional
        Type of minor helix, can be: 'alpha', 'pi', '3-10',
        'PPI', 'PPII', 'collagen'.
    orientation : int, optional
        Parallel (1) or anti-parallel (-1).
    phi_c_alpha : float, optional
        Rotation of the minor helix relative to the super-
        helical axis.
    minor_repeat : float, optional
        Hydrophobic repeat of the `Helix`.

    Attributes
    ----------
    num_monomers : int
        Number of amino acids in the `Helix`. If `None, an
    major_pitch : float, optional
        Pitch of the super helix.
    major_radius : float, optional
        Radius of the super helix.
    major_handedness : str, optional
        Handedness of the super helix.
    minor_helix_type : str, optional
        Type of minor helix, can be: 'alpha', 'pi', '3-10',
        'PPI', 'PPII', 'collagen'.
    orientation : int, optional
        Parallel (1) or anti-parallel (-1).
    phi_c_alpha : float, optional
        Rotation of the minor helix relative to the super-
        helical axis.
    minor_repeat : float, optional
        Hydrophobic repeat of the `Helix`,
        appropriate number of residues are added.
    minor_rise_per_residue : float
        Distance that each residue travels along the z-axis.
    minor_handedness : str
        The handedness of the helix can be 'r' or 'l'.
    helix_start : 3D Vector (tuple or list or numpy.array)
        The coordinate of the start of the helix primitive.
    helix_end : 3D Vector (tuple or list or numpy.array)
        The coordinate of the end of the helix primitive.
    """

    def __init__(self, aa=10, major_pitch=225.8, major_radius=5.07,
                 major_handedness='l', minor_helix_type='alpha', orientation=1,
                 phi_c_alpha=0.0, minor_repeat=None):
        super(HelicalHelix, self).__init__()
        # Major helix properties
        self.num_monomers = aa
        self.major_pitch = major_pitch
        self.major_radius = major_radius
        self.major_handedness = major_handedness
        # Minor helix properties
        self.minor_helix_type = minor_helix_type
        self.minor_rise_per_residue = _helix_parameters[self.minor_helix_type][1]
        self.minor_handedness = _helix_handedness[self.minor_helix_type]
        self.orientation = orientation
        if minor_repeat == 0:
            minor_repeat = None
        self.minor_repeat = minor_repeat
        self.phi_c_alpha = phi_c_alpha
        # major helix start and end.
        self.helix_start = numpy.array([0.0, 0.0, 0.0])
        self.helix_end = self.helix_length * numpy.array([0.0, 0.0, 1.0])
        self.build()

    def __repr__(self):
        if len(self.sequence) > 15:
            seq = self.sequence[:12] + '...'
        else:
            seq = self.sequence
        return '<HelicalHelix containing {} {}. Sequence: {}>'.format(
            len(self._monomers), 'Residue' if len(self._monomers) == 1 else 'Residues', seq)

    @classmethod
    def from_start_and_end(cls, start, end, aa=None, major_pitch=225.8,
                           major_radius=5.07, major_handedness='l',
                           minor_helix_type='alpha', orientation=1,
                           phi_c_alpha=0.0, minor_repeat=None):
        """Creates a `HelicalHelix` between a `start` and `end` point."""
        start = numpy.array(start)
        end = numpy.array(end)
        if aa is None:
            minor_rise_per_residue = _helix_parameters[minor_helix_type][1]
            aa = int((numpy.linalg.norm(end - start) /
                      minor_rise_per_residue) + 1)
        instance = cls(
            aa=aa, major_pitch=major_pitch, major_radius=major_radius,
            major_handedness=major_handedness,
            minor_helix_type=minor_helix_type, orientation=orientation,
            phi_c_alpha=phi_c_alpha, minor_repeat=minor_repeat)
        instance.move_to(start=start, end=end)
        return instance

    @property
    def helix_length(self):
        """ Length of major axis (i.e. height of cylinder around which helicalhelix curves)."""
        return self.num_monomers * self.major_rise_per_monomer

    @property
    def major_axis(self):
        """Axis of the super helix."""
        return Axis(start=self.helix_start, end=self.helix_end)

    @property
    def curve(self):
        """Curve of the super helix."""
        return HelicalCurve.pitch_and_radius(
            self.major_pitch, self.major_radius,
            handedness=self.major_handedness)

    @property
    def curve_primitive(self):
        """`Primitive` of the super-helical curve."""
        curve = self.curve
        curve.axis_start = self.helix_start
        curve.axis_end = self.helix_end
        coords = curve.get_coords(
            n_points=(self.num_monomers + 1), spacing=self.minor_rise_per_residue)
        if self.orientation == -1:
            coords.reverse()
        return Primitive.from_coordinates(coords)

    @property
    def major_rise_per_monomer(self):
        """Rise along super-helical axis per monomer."""
        return numpy.cos(numpy.deg2rad(self.curve.alpha)) * self.minor_rise_per_residue

    def minor_residues_per_turn(self, minor_repeat=None):
        """Calculates the number of residues per turn of the minor helix.

        Parameters
        ----------
        minor_repeat : float, optional
            Hydrophobic repeat of the minor helix.

        Returns
        -------
        minor_rpt : float
            Residues per turn of the minor helix.
        """
        if minor_repeat is None:
            minor_rpt = _helix_parameters[self.minor_helix_type][0]
        else:
            # precession angle in radians
            precession = self.curve.t_from_arc_length(
                minor_repeat * self.minor_rise_per_residue)
            if self.orientation == -1:
                precession = -precession
            if self.major_handedness != self.minor_handedness:
                precession = -precession
            minor_rpt = ((minor_repeat * numpy.pi * 2) /
                         ((2 * numpy.pi) + precession))
        return minor_rpt

    def build(self):
        """Builds the `HelicalHelix`."""
        helical_helix = Polypeptide()
        primitive_coords = self.curve_primitive.coordinates
        helices = [Helix.from_start_and_end(start=primitive_coords[i],
                                            end=primitive_coords[i + 1],
                                            helix_type=self.minor_helix_type,
                                            aa=1)
                   for i in range(len(primitive_coords) - 1)]
        residues_per_turn = self.minor_residues_per_turn(
            minor_repeat=self.minor_repeat)
        if residues_per_turn == 0:
            residues_per_turn = _helix_parameters[self.minor_helix_type][0]
        if self.minor_handedness == 'l':
            residues_per_turn *= -1
        # initial phi_c_alpha value calculated using the first Helix in helices.
        if self.orientation != -1:
            initial_angle = dihedral(numpy.array([0, 0, 0]),
                                     primitive_coords[0],
                                     primitive_coords[1],
                                     helices[0][0]['CA'])
        else:
            initial_angle = dihedral(
                numpy.array([0, 0, primitive_coords[0][2]]),
                primitive_coords[0],
                numpy.array([primitive_coords[0][0],
                             primitive_coords[0][1], primitive_coords[1][2]]),
                helices[0][0]['CA'])
        # angle required to achieve desired phi_c_alpha value of self.phi_c_alpha.
        addition_angle = self.phi_c_alpha - initial_angle
        for i, h in enumerate(helices):
            angle = (i * (360.0 / residues_per_turn)) + addition_angle
            h.rotate(angle=angle, axis=h.axis.unit_tangent,
                     point=h.helix_start)
            helical_helix.extend(h)
        helical_helix.relabel_all()
        self._monomers = helical_helix._monomers[:]
        for monomer in self._monomers:
            monomer.parent = self
        return

    def get_orient_angle(self, reference_point=numpy.array([0, 0, 0]),
                         monomer_index=0, res_label='CA', radians=False):
        """ Angle between reference_point and self[monomer_index][res_label].

        Notes
        -----
        Angle is calculated using the dihedral angle, with the 
        second and third points coming from the curve_primitive.

        Parameters
        ----------
        reference_point : list, tuple or numpy.array of length 3.
        monomer_index : int
            Index of the Residue to centre.
        res_label : str
            Atom name for centred atom, e.g. "CA" or "OE1".
        radians : bool
            If True, then desired_angle is in radians instead of degrees.
        """
        if (monomer_index < len(self)) and monomer_index != -1:
            adjacent_index = monomer_index + 1
        elif (monomer_index == len(self)) or monomer_index == -1:
            adjacent_index = monomer_index - 1
        else:
            raise ValueError(
                "centred_index ({0}) cannot be greater than the "
                "length of the polymer ({1})".format(
                    monomer_index, len(self)))
        angle = dihedral(reference_point,
                         self.curve_primitive[monomer_index]['CA'],
                         self.curve_primitive[adjacent_index]['CA'],
                         self[monomer_index][res_label])
        if radians:
            angle = numpy.deg2rad(angle)
        return angle

    def translate(self, vector):
        super(HelicalHelix, self).translate(vector=vector)
        self.helix_start += vector
        self.helix_end += vector
        return

    def rotate(self, angle, axis, point=None, radians=False):
        super(HelicalHelix, self).rotate(angle=angle,
                                         axis=axis, point=point, radians=radians)
        # modify helix_start and helix_end accordingly.
        q = Quaternion.angle_and_axis(angle=angle, axis=axis, radians=radians)
        self.helix_start = q.rotate_vector(v=self.helix_start, point=point)
        self.helix_end = q.rotate_vector(v=self.helix_end, point=point)
        return

    def move_to(self, start, end):
        """Moves the `HelicalHelix` to lie on the vector `start` and `end`.

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

    def rotate_monomers(self, angle, radians=False):
        """ Rotates each Residue in the Polypeptide.

        Notes
        -----
        Each monomer is rotated about the axis formed between its
        corresponding primitive `PseudoAtom` and that of the 
        subsequent `Monomer`.

        Parameters
        ----------
        angle : float
            Angle by which to rotate each monomer.
        radians : bool
            Indicates whether angle is in radians or degrees.
        """
        if radians:
            angle = numpy.rad2deg(angle)
        for i in range(len(self.primitive) - 1):
            axis = self.primitive[i + 1]['CA'] - self.primitive[i]['CA']
            point = self.primitive[i]['CA']._vector
            self[i].rotate(angle=angle, axis=axis, point=point)
        return


__author__ = 'Jack W. Heal, Christopher W. Wood'
