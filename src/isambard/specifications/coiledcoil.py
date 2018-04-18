"""Contains code for modelling coiled coils and collagens."""

from ampal import Assembly
import numpy

from .helix import HelicalHelix, _helix_parameters

basis_set_parameters = {
    2: {'name': 'CC-Di', 'pitch': 225.8, 'radius': 5.07, 'interface_angle': 283.56,
        'sequence': 'EIAALKQEIAALKKENAALKWEIAALKQ'},
    3: {'name': 'CC-Tri', 'pitch': 194.0, 'radius': 6.34, 'interface_angle': 277.12,
        'sequence': 'EIAAIKQEIAAIKKEIAAIKWEIAAIKQ'},
    4: {'name': 'CC-Tet', 'pitch': 213.2, 'radius': 6.81, 'interface_angle': 279.20,
        'sequence': 'ELAAIKQELAAIKKELAAIKWELAAIKQ'},
    5: {'name': 'CC-Pent', 'pitch': 182.8, 'radius': 8.62, 'interface_angle': 271.58,
        'sequence': 'KIEQILQKIEKILQKIEWILQKIEQILQ'},
    6: {'name': 'CC-Hex', 'pitch': 228.4, 'radius': 9.13, 'interface_angle': 273.28,
        'sequence': 'ELKAIAQELKAIAKELKAIAWELKAIAQ'},
    7: {'name': 'CC-Hept', 'pitch': 328.6, 'radius': 9.80, 'interface_angle': 272.24,
        'sequence': 'EIAQALKEIAKALKEIAWALKEIAQALK'},
}


class CoiledCoil(Assembly):
    """Models a coiled-coil protein.

    Notes
    -----
    Instantiating this class using just an oligomeric state is used
    to create simple reference models. To build more complex models
    use the `from_parameters` classmethod.

    Parameters
    ----------
    n : int
        The oligomeric state of the model to be built.
    auto_build : bool, optional
        If `True`, the model will be built as part of instantiation.

    Attributes
    ----------
    aas : [int]
        Number of amino acids in each minor helix.
    basis_set_sequences : [str]
        Reference sequences for the oligomeric state that has been
        selected, taken from the basis set of coiled coils.
    major_radii : [float]
        Radii of the minor helices relative to the super-helical
        axis.
    major_pitches : [float]
        Pitch values of the minor helices relative to the super-helical
        axis.
    phi_c_alphas :
        Relative rotation values of the minor helices relative to
        the super-helical axis.
    minor_helix_types : [str]
        Helix types of the minor helices. Can be: 'alpha', 'pi', '3-10',
        'PPI', 'PP2', 'collagen'.
    major_handedness : str
        Handedness of the super helix.
    orientations :
        Orientation of helices relative to the super-helical axis. 1
        is parallel, -1 is anti-parallel.
    minor_repeats : [float]
        Hydrophobic repeats of the minor helices.
    rotational_offsets :
        Rotation of the minor helices relative to the super-helical
        axis.
    z_shifts : [float]
        Translation of the minor helices along the super-helical axis.
    oligomeric_state : int
        Oligomeric state of the coiled coil.
    """

    def __init__(self, n, auto_build=True):
        super(CoiledCoil, self).__init__()
        # parameters for each polypeptide
        # basis set parameters if known, otherwise educated guesses.
        if n in basis_set_parameters.keys():
            parameters = basis_set_parameters[n]
            radius = parameters['radius']
        else:
            # calculate radius based on extrapolated straight-line fit
            # of n Vs radius for basis_set_parameters
            radius = (n * 0.966) + 3.279
            # other default values just copied from largest oligomer
            # in basis_set_parameters.
            parameters = basis_set_parameters[max(basis_set_parameters.keys())]
        self.major_radii = [radius] * n
        self.major_pitches = [parameters['pitch']] * n
        self.basis_set_sequences = [parameters['sequence']] * n
        self.aas = [len(parameters['sequence'])] * n
        self.phi_c_alphas = [parameters['interface_angle']] * n
        # alpha-helical barrel with heptad repeat as default.
        self.major_handedness = ['l'] * n
        self.minor_helix_types = ['alpha'] * n
        self.orientations = [1] * n
        self.minor_repeats = [3.5] * n
        # parameters for the arrangement of each polypeptide
        # (evenly distributed, no z-displacement).
        self.rotational_offsets = [((i * 360.0) / n) for i in range(n)]
        self.z_shifts = [0.0] * n
        # parameters for the whole assembly
        self.oligomeric_state = n
        if auto_build:
            self.build()

    @classmethod
    def from_polymers(cls, polymers):
        """Creates a `CoiledCoil` from a list of `HelicalHelices`.

        Parameters
        ----------
        polymers : [HelicalHelix]
            List of `HelicalHelices`.
        """
        n = len(polymers)
        instance = cls(n=n, auto_build=False)
        instance.major_radii = [x.major_radius for x in polymers]
        instance.major_pitches = [x.major_pitch for x in polymers]
        instance.major_handedness = [x.major_handedness for x in polymers]
        instance.aas = [x.num_monomers for x in polymers]
        instance.minor_helix_types = [x.minor_helix_type for x in polymers]
        instance.orientations = [x.orientation for x in polymers]
        instance.phi_c_alphas = [x.phi_c_alpha for x in polymers]
        instance.minor_repeats = [x.minor_repeat for x in polymers]
        instance.build()
        return instance

    @classmethod
    def from_parameters(cls, n, aa=28, major_radius=None, major_pitch=None,
                        phi_c_alpha=26.42, minor_helix_type='alpha',
                        auto_build=True):
        """Creates a `CoiledCoil` from defined super-helical parameters.

        Parameters
        ----------
        n : int
            Oligomeric state
        aa : int, optional
            Number of amino acids per minor helix.
        major_radius : float, optional
            Radius of super helix.
        major_pitch : float, optional
            Pitch of super helix.
        phi_c_alpha : float, optional
            Rotation of minor helices relative to the super-helical
            axis.
        minor_helix_type : float, optional
            Helix type of minor helices. Can be: 'alpha', 'pi', '3-10',
            'PPI', 'PP2', 'collagen'.
        auto_build : bool, optional
            If `True`, the model will be built as part of instantiation.
        """
        instance = cls(n=n, auto_build=False)
        instance.aas = [aa] * n
        instance.phi_c_alphas = [phi_c_alpha] * n
        instance.minor_helix_types = [minor_helix_type] * n
        if major_pitch is not None:
            instance.major_pitches = [major_pitch] * n
        if major_radius is not None:
            instance.major_radii = [major_radius] * n
        if auto_build:
            instance.build()
        return instance

    @classmethod
    def tropocollagen(
            cls, aa=28, major_radius=5.0, major_pitch=85.0, auto_build=True):
        """Creates a model of a collagen triple helix.

        Parameters
        ----------
        aa : int, optional
            Number of amino acids per minor helix.
        major_radius : float, optional
            Radius of super helix.
        major_pitch : float, optional
            Pitch of super helix.
        auto_build : bool, optional
            If `True`, the model will be built as part of instantiation.
        """
        instance = cls.from_parameters(
            n=3, aa=aa, major_radius=major_radius, major_pitch=major_pitch,
            phi_c_alpha=0.0, minor_helix_type='collagen', auto_build=False)
        instance.major_handedness = ['r'] * 3
        # default z-shifts taken from rise_per_residue of collagen helix
        rpr_collagen = _helix_parameters['collagen'][1]
        instance.z_shifts = [-rpr_collagen * 2, -rpr_collagen, 0.0]
        instance.minor_repeats = [None] * 3
        if auto_build:
            instance.build()
        return instance

    def build(self):
        """Builds a model of a coiled coil protein using input parameters."""
        monomers = [HelicalHelix(major_pitch=self.major_pitches[i],
                                 major_radius=self.major_radii[i],
                                 major_handedness=self.major_handedness[i],
                                 aa=self.aas[i],
                                 minor_helix_type=self.minor_helix_types[i],
                                 orientation=self.orientations[i],
                                 phi_c_alpha=self.phi_c_alphas[i],
                                 minor_repeat=self.minor_repeats[i],
                                 )
                    for i in range(self.oligomeric_state)]
        axis_unit_vector = numpy.array([0, 0, 1])
        for i, m in enumerate(monomers):
            m.rotate(angle=self.rotational_offsets[i], axis=axis_unit_vector)
            m.translate(axis_unit_vector * self.z_shifts[i])
        self._molecules = monomers[:]
        self.relabel_all()
        for m in self._molecules:
            m.ampal_parent = self
        return


__author__ = 'Andrew R. Thomson, Christopher W. Wood, Jack W. Heal'
