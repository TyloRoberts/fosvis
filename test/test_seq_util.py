from fosvis import seq_util

import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import filecmp
import os
import pandas as pd
from pandas.testing import assert_frame_equal
import difflib


# Done
class test_remove_too_small_contigs(unittest.TestCase):

    def test_remove_too_small_contigs(self):
        input_file = 'test/data_for_testing/Fosmid_Size_Selection_Tests/test_remove_too_small_contigs_input.fasta'
        result = seq_util.remove_too_small_contigs(input_file, 100)

        seq3_len_100 = SeqRecord(
            Seq("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
            id="seq3_len_100",
            name="seq3_len_100",
            description="seq3_len_100")
        seq4_len_101 = SeqRecord(
            Seq("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
            id="seq4_len_101",
            name="seq4_len_101",
            description="seq4_len_101")
        seq5_len_200 = SeqRecord(
            Seq("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
            id="seq5_len_200",
            name="seq5_len_200",
            description="seq5_len_200")
        expected_result = [seq3_len_100, seq4_len_101, seq5_len_200]

        self.assertEqual(len(result), len(expected_result))
        self.assertTrue(result[0].seq == expected_result[0].seq)
        self.assertTrue(result[1].seq == expected_result[1].seq)
        self.assertTrue(result[2].seq == expected_result[2].seq)
        self.assertTrue(result[0].id == expected_result[0].id)
        self.assertTrue(result[1].id == expected_result[1].id)
        self.assertTrue(result[2].id == expected_result[2].id)

# Done
# class test_write_seqs_to_file(unittest.TestCase):
#
#     def test_write_seqs_to_file(self):
#         file_to_write = 'test/data_for_testing/Fosmid_Size_Selection_Tests/write_seqs_to_file_actual_output.fasta'
#         seq1 = SeqRecord(
#             Seq("AAGGTTCC"),
#             id="seq1",
#             name="seq1",
#             description="seq1")
#         seq2 = SeqRecord(
#             Seq("GGAACCTT"),
#             id="seq2",
#             name="seq2",
#             description="seq2")
#         seqs_to_write = [seq1,seq2]
# 
#         seq_util.write_seqs_to_file(seqs_to_write, file_to_write)
#
#         expected_output_file = 'test/data_for_testing/Fosmid_Size_Selection_Tests/test_write_seqs_to_file_expected_output.fasta'
#
#         self.assertTrue(filecmp.cmp(expected_output_file, file_to_write, shallow=False))
#
#         os.remove(file_to_write)

# Done
class test_get_karyotype_data(unittest.TestCase):

    def test_get_karyotype_data(self):
        seq1_len_8 = SeqRecord(
            Seq("AAGGTTCC"),
            id="seq1_len_8",
            name="seq1_len_8",
            description="seq1_len_8")
        seq2_len_18 = SeqRecord(
            Seq("GGAACCTTGGAACCTT"),
            id="seq1_len_8",
            name="seq1_len_8",
            description="seq1_len_8")
        seqs = [seq1_len_8,seq2_len_18]
        result = seq_util.get_karyotype_data(seqs)

        expected_result_data = {'chr_prefix':['chr', 'chr'], '-prefix':['-', '-'], 'variable_name':['seq1_len_8', 'seq1_len_8'], 'diagram_label':['1', '2'], 'start':[1, 1], 'end':[9, 17], 'color':['rgb(120,120,120,0.4)', 'rgb(120,120,120,0.4)']}
        expected_result_df = pd.DataFrame(expected_result_data)

        self.assertEqual(assert_frame_equal(result, expected_result_df, check_dtype=False), None)

# Done
class test_gc_interval(unittest.TestCase):

    def test_gc_interval(self):

        interval_10_input_file = 'test/data_for_testing/gc_interval_tests/gc_content_test_interval_10_input.fasta'
        interval_3_input_file = 'test/data_for_testing/gc_interval_tests/gc_content_test_interval_3_input.fasta'

        interval_10_result = seq_util.gc_interval(interval_10_input_file, 10)
        interval_3_result = seq_util.gc_interval(interval_3_input_file, 3)

        expected_interval_10_result = {'contig':['interval_10_mix', 'interval_10_mix', 'interval_10_mix', 'interval_10_mix',
                                                 'interval_10_all_GC', 'interval_10_all_GC', 'interval_10_not_divis_by_10'],
                                       'interval_start':[1, 11, 21, 31, 1, 11, 1], 'interval_end':[11,21,31,41, 11, 21, 5], 'gc_content':[100,50,50,0,100,100,50]}
        expected_interval_10_result_df = pd.DataFrame(expected_interval_10_result)
        expected_interval_3_result = {'contig':['interval_3_mix', 'interval_3_mix', 'interval_3_mix', 'interval_3_mix',
                                                 'interval_3_not_divis_by_3', 'interval_3_not_divis_by_3', 'interval_3_not_divis_by_3', 'interval_3_not_divis_by_3'],
                                       'interval_start':[1,4,7,10, 1, 4, 7, 10],
                                       'interval_end':[4, 7, 10, 13, 4, 7, 10, 11],
                                       'gc_content':[((2/3)*100), 0, 100, (1/3)*100, 100,100,0,100]}

        expected_interval_3_result_df = pd.DataFrame(expected_interval_3_result)

        self.assertEqual(assert_frame_equal(interval_10_result, expected_interval_10_result_df, check_dtype=False), None)
        self.assertEqual(assert_frame_equal(interval_3_result, expected_interval_3_result_df, check_dtype=False), None)


if __name__ == "__main__":
    unittest.main()
