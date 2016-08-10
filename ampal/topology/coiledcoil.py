import numpy

from ampal.secondary_structure.helix import HelicalHelix, _helix_parameters
from ampal.protein import Assembly

basis_set_parameters = {
    2: {'name': 'CC-Di', 'pitch': 225.8, 'radius': 5.07, 'interface_angle': 26.42,
        'sequence': 'EIAALKQEIAALKKENAALKWEIAALKQ'},
    3: {'name': 'CC-Tri', 'pitch': 194.0, 'radius': 6.34, 'interface_angle': 19.98,
        'sequence': 'EIAAIKQEIAAIKKEIAAIKWEIAAIKQ'},
    4: {'name': 'CC-Tet', 'pitch': 213.2, 'radius': 6.81, 'interface_angle': 22.06,
        'sequence': 'ELAAIKQELAAIKKELAAIKWELAAIKQ'},
    5: {'name': 'CC-Pent', 'pitch': 182.8, 'radius': 8.62, 'interface_angle': 14.44,
        'sequence': 'KIEQILQKIEKILQKIEWILQKIEQILQ'},
    6: {'name': 'CC-Hex', 'pitch': 228.4, 'radius': 9.13, 'interface_angle': 16.14,
        'sequence': 'ELKAIAQELKAIAKELKAIAWELKAIAQ'},
    7: {'name': 'CC-Hept', 'pitch': 328.6, 'radius': 9.80, 'interface_angle': 15.10,
        'sequence': 'EIAQALKEIAKALKEIAWALKEIAQALK'},
}


class CoiledCoil(Assembly):
    def __init__(self, n, auto_build=True):
        super(CoiledCoil, self).__init__()
        # parameters for each polypeptide
        # basis set parameters if known, otherwise educated guesses.
        if n in basis_set_parameters.keys():
            parameters = basis_set_parameters[n]
            radius = parameters['radius']
        else:
            # calculate radius based on extrapolated straight-line fit of n Vs radius for basis_set_parameters
            radius = (n * 0.966) + 3.279
            # other default values just copied from largest oligomer in basis_set_parameters.
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
        # parameters for the arrangement of each polypeptide (evenly distributed, no z-displacement).
        self.rotational_offsets = [((i * 360.0) / n) for i in range(n)]
        self.z_shifts = [0.0] * n
        # parameters for the whole assembly
        self.oligomeric_state = n
        if auto_build:
            self.build()

    @classmethod
    def from_polymers(cls, polymers):
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
    def from_parameters(cls, n, aa=28, major_radius=None, major_pitch=None, phi_c_alpha=26.42,
                        minor_helix_type='alpha', auto_build=True):
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
    def tropocollagen(cls, aa=28, major_radius=5.0, major_pitch=85.0, auto_build=True):
        instance = cls.from_parameters(n=3, aa=aa, major_radius=major_radius, major_pitch=major_pitch,
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


