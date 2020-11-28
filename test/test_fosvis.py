import sys, os
testdir = os.path.dirname(__file__)
srcdir = '../fosvis'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))

import fosvis
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import filecmp
import pandas as pd
from pandas.testing import assert_frame_equal
import difflib


"""
To run all: python3 -m unittest
To run a specific test: python3 -m unittest test.test_fosvis
"""

# Done
class test_setup_project(unittest.TestCase):

    def test_setup_project_directory_already_exists(self):
        """
        Tests that when trying to setup a project with the same directory that
        already exists it will throw a SystemExit error.
        """
        test_dir_path = 'test/data_for_testing/testing_temporary_files/setup_project_test_dir'
        os.mkdir(test_dir_path)
        with self.assertRaises(SystemExit):
            fosvis.setup_project_dirs(test_dir_path, 'blastn')
        os.rmdir(test_dir_path)


class test_create_circos_data(unittest.TestCase):

    def test_(self):
        pass


class test_write_paramaters_log(unittest.TestCase):

    def test_write_paramaters_log(self):
        pass


if __name__ == "__main__":
    unittest.main()
