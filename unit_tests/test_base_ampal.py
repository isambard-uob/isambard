from collections import Counter
import os
import unittest

import isambard_dev as isambard


class PDBWriterTestCase(unittest.TestCase):
    """Tests for isambard.tools.tool_geometry.dihedral"""

    def pdb_check(self, pdb_file_path):
        with open(pdb_file_path, 'r') as inf:
            pdb_file = inf.read()
        ampal = isambard.ampal.convert_pdb_to_ampal(pdb_file_path)

        # Compare the number of lines in the output ampal pdb with the original
        pdb_lines = [x for x in pdb_file.splitlines() if (x.startswith('ATOM') or x.startswith('HETATM'))]
        ampal_pdb_lines = [x for x in ampal.make_pdb(
                ligands=True, alt_states=True).splitlines() if (x.startswith('ATOM') or x.startswith('HETATM'))]
        self.assertEqual(len(pdb_lines), len(ampal_pdb_lines))

        # Compare the atomic composition
        pdb_atomic = Counter([x[-4:-2].strip() for x in pdb_lines])
        ampal_atomic = Counter([x[-4:-2].strip() for x in ampal_pdb_lines])
        self.assertEqual(pdb_atomic, ampal_atomic)

        # Compare the residue composition
        pdb_atomic = Counter([x[17:20].strip() for x in pdb_lines])
        ampal_atomic = Counter([x[17:20].strip() for x in ampal_pdb_lines])
        self.assertEqual(pdb_atomic, ampal_atomic)
        return

    def test_3qy1(self):
        """Simple test for reading in protein and ligands."""
        test_file_path = os.path.join('unit_tests', 'testing_files', '3qy1.pdb')
        self.pdb_check(test_file_path)

    def test_2ht0(self):
        """More complex test for reading in protein, DNA and ligands."""
        test_file_path = os.path.join('unit_tests', 'testing_files', '2ht0.pdb')
        self.pdb_check(test_file_path)

    def test_1ek9(self):
        """Very complex pdb file of large protein channel."""
        test_file_path = os.path.join('unit_tests', 'testing_files', '1ek9.pdb')
        self.pdb_check(test_file_path)

if __name__ == '__main__':
    unittest.main()
