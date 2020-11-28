from fosvis import diagram_and_legend

import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import filecmp
import os
import pandas as pd
from pandas.testing import assert_frame_equal
import difflib


class test_create_circos_conf_file(unittest.TestCase):

    def test_create_circos_conf_file_no_gc(self):
        no_gc_actual_ouput_path = 'test/data_for_testing/test_diagram_and_legend/test_create_circos_conf_file/circos_actual_no_gc.conf'
        diagram_and_legend.create_circos_conf_file(no_gc_actual_ouput_path, 'karyotype_path.txt', 'links_path.txt', 'ORF_path.txt', 'ORF_reverse_path.txt')
        expected_output_file_no_gc = 'test/data_for_testing/test_diagram_and_legend/test_create_circos_conf_file/circos_expected_no_gc.conf'
        self.assertTrue(filecmp.cmp(expected_output_file_no_gc, no_gc_actual_ouput_path, shallow=False))
        os.remove(no_gc_actual_ouput_path)


    def test_create_circos_conf_file_with_gc(self):
        with_gc_actual_ouput_path = 'test/data_for_testing/test_diagram_and_legend/test_create_circos_conf_file/circos_actual_with_gc.conf'
        diagram_and_legend.create_circos_conf_file(with_gc_actual_ouput_path, 'karyotype_path.txt', 'links_path.txt', 'ORF_path.txt', 'ORF_reverse_path.txt', gc_data='gc_content.txt')
        expected_output_file_with_gc = 'test/data_for_testing/test_diagram_and_legend/test_create_circos_conf_file/circos_expected_with_gc.conf'
        self.assertTrue(filecmp.cmp(expected_output_file_with_gc, with_gc_actual_ouput_path, shallow=False))
        os.remove(with_gc_actual_ouput_path)



class test_make_diagram(unittest.TestCase):

    def test_(self):
        pass


class test_prtn_domain_legend(unittest.TestCase):

    def test_(self):
        pass


class test_draw_prtn_domain_legend(unittest.TestCase):

    def test_(self):
        pass


class test_karyotype_legend(unittest.TestCase):

    def test_(self):
        pass


if __name__ == "__main__":
    unittest.main()
