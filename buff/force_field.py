import json
import os

from settings import global_settings

force_fields = {}
for ff in os.listdir(os.path.join(global_settings['package_path'], 'buff', 'force_fields')):
    ffs = ff.split('.')
    if ffs[-1] == 'json':
        force_fields[ffs[0]] = os.path.join(global_settings['package_path'], 'buff', 'force_fields', ff)


class BuffForceField(dict):
    def __init__(self, force_field='standard'):
        with open(force_fields[force_field], 'r') as inf:
            in_d = json.loads(inf.read())
        super().__init__(in_d)
        self.force_field = force_field

    def __repr__(self):
        return "<BUFF Force Field Object: {}>".format(self.force_field)

    @property
    def max_radius_and_npnp(self):
        return self.find_max_rad_npnp()

    @property
    def distance_cutoff(self):
        rad, npnp = self.find_max_rad_npnp()
        return (rad * 2) + npnp

    def find_max_rad_npnp(self):
        max_rad = 0
        max_npnp = 0
        for res, atoms in self.items():
            if res != 'KEY':
                for atom, ff_params in self[res].items():
                    if max_rad < ff_params[1]:
                        max_rad = ff_params[1]
                    if max_npnp < ff_params[4]:
                        max_npnp = ff_params[4]
        return max_rad, max_npnp

global_settings[u'buff'][u'force_field'] = BuffForceField(
    force_field=global_settings[u'buff'][u'default_force_field'])
