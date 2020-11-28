from fosvis import orf

import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import filecmp
import os
import pandas as pd
from pandas.testing import assert_frame_equal
import difflib

# Done
class test_run_prodigal(unittest.TestCase):
    """
    Checks that prodigal creates both output files
    Does not check that contents of the file
    """

    def test_run_prodigal(self):
        sample_seq_path = './test/data_for_testing/test_orf_data/test_orf_sample_seq.fasta'
        output_path = './test/data_for_testing/test_orf_data/temp_outputs/orf_test_output.txt'
        prtn_output_path = './test/data_for_testing/test_orf_data/temp_outputs/orf_test_prtn_output.faa'
        orf.run_prodigal(sample_seq_path, output_path, prtn_output_path)

        # Check both outputs exist and are non-empty
        self.assertTrue(os.path.exists(output_path) and os.path.getsize(output_path) > 0)
        self.assertTrue(os.path.exists(prtn_output_path) and os.path.getsize(prtn_output_path) > 0)

        # remove the files
        os.remove(output_path)
        os.remove(prtn_output_path)


# Done
class test_parse_prodigal_sco_output(unittest.TestCase):

    def test_parse_prodigal_sco_output(self):
        # Parse a sample file
        prodigal_sample_output_path = './test/data_for_testing/test_orf_data/test_parse_prodigal_sco_output_input.txt'
        output = orf.parse_prodigal_sco_output(prodigal_sample_output_path)
        forward_strand_genes_actual_df = output[0]
        reverse_strand_genes_actual_df = output[1]

        #  Create the expected output
        forward_strand_genes_expected_df = pd.DataFrame({'contig' : ['contig_1', 'contig_1', 'contig_1', 'contig_1', 'contig_2', 'contig_2', 'contig_2', 'contig_2', 'contig_2'],
                                                         'orf_start' : [1,759,1688,2361,18443,19870,20100,22725,23899],
                                                         'orf_end' : [547,1686,2345,3192,19451,20077,20574,23064,25633]})


        reverse_strand_genes_expected_df = pd.DataFrame({'contig' : ['contig_1', 'contig_1', 'contig_1', 'contig_1', 'contig_1', 'contig_1', 'contig_2', 'contig_2', 'contig_2', 'contig_2', 'contig_2', 'contig_2', 'contig_2', 'contig_2', 'contig_2', 'contig_2'],
                                                         'orf_start' : [3228,4094,5473,7205,7623,8844,95,1768,2273,11474,12047,13858,14829,15125,17691,17807],
                                                         'orf_end' : [4095,5477,7108,7532,8835,9141,1145,2245,2999,11777,12371,14521,15015,16082,17811,18038]})


        # Check that actual and expected have same shape
        self.assertEqual(forward_strand_genes_expected_df.shape, forward_strand_genes_actual_df.shape)
        self.assertEqual(reverse_strand_genes_expected_df.shape, reverse_strand_genes_actual_df.shape)

        # Check contents is equal
        self.assertEqual(assert_frame_equal(forward_strand_genes_expected_df, forward_strand_genes_actual_df), None)
        self.assertEqual(assert_frame_equal(reverse_strand_genes_expected_df, reverse_strand_genes_actual_df), None)




if __name__ == "__main__":
    unittest.main()
