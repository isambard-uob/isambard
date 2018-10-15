
import pathlib
import unittest

import isambard
import ampal
import isambard.evaluation.packing_density as packing_density

TEST_FILES_PATH = pathlib.Path(__file__).parent / 'testing_files'

class PackingDensityTestCase(unittest.TestCase):
    """Tests for isambard.evaluation.packing_density"""

    def setUp(self):
        self.ht0 = ampal.load_pdb(str(TEST_FILES_PATH / '2ht0.pdb'))
        self.ubq = ampal.load_pdb(str(TEST_FILES_PATH / '1ubq.pdb'))
        self.pdbs = [self.ht0, self.ubq]

    def check_packing_density(self):
        for pdb in self.pdbs:
            # Creates dictionary of packing density values from ISAMBARD function
            atoms_list_func = list(packing_density.calculate_packing_density(pdb))
            pack_dens_func_val = {atom.id: atom.tags['packing density'] for
                                  atom in atoms_list_func}
            pack_dens_func_val = sorted(pack_dens_func_val.items())

            # Creates dictionary of packing density values from is_within AMPAL
            # function
            atoms_list_calc = [atom for atom in pdb.get_atoms()
                               if atom.element != 'H']
            pack_dens_calc_val = {}
            for atom in atoms_list_calc:
                pack_dens_calc_val[atom.id] = (len(pdb.is_within(7, atom)) - 1)
            pack_dens_calc_val = sorted(pack_dens_calc_val.items())

            self.assertEqual(pack_dens_func_val, pack_dens_calc_val)
