from fosvis import protein_domains

import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import filecmp
import os
import pandas as pd
from pandas.testing import assert_frame_equal
import difflib


class test_run_and_parse_hmmscan(unittest.TestCase):

    def test_run_and_parse_hmmscan(self):
        hmm_db = '/mnt/nfs/sharknado/Sandbox/Tylo/databases/hmmer/hmmscan/Pfam-A.hmm'

        # Run run_hmmscan() on sample orfs
        sample_orfs = './test/data_for_testing/testing_protein_domains/sample_orfs_and_hmmscan_df/prodigal_prtn_seq_output.faa'
        output_dir = './test/data_for_testing/testing_protein_domains/temp_hmmscan_output'
        protein_domains.run_hmmscan(sample_orfs, hmm_db, output_dir, 0.01, 'hmmscan_testing')

        # check output file exits and is not empty
        tblout_path = output_dir + '/hmmscan_testing_tblout.txt'
        domtblout_path = output_dir + '/hmmscan_testing_domtblout.txt'
        out_file_path = output_dir + '/hmmscan_testing_out.txt'
        self.assertTrue(os.path.exists(tblout_path) and os.path.getsize(tblout_path) > 0)
        self.assertTrue(os.path.exists(domtblout_path) and os.path.getsize(domtblout_path) > 0)
        self.assertTrue(os.path.exists(out_file_path) and os.path.getsize(out_file_path) > 0)

        # Parse the output files using parse_hmmscan_domtblout()
        actual_output_df = protein_domains.parse_hmmscan_domtblout(domtblout_path)

        # Check df
        self.assertEqual(5, actual_output_df.shape[1])
        self.assertEqual(['target', 'orf_id', 'orf_start', 'orf_stop', 'target_description'], list(actual_output_df.columns.values))

        # Delete all files
        os.remove(tblout_path)
        os.remove(domtblout_path)
        os.remove(out_file_path)


class test_convert_hmmscan_domains_from_protein_to_nuc_pos(unittest.TestCase):

    def test_convert_hmmscan_domains_from_protein_to_nuc_pos(self):
        sample_hmmscan_df_path = './test/data_for_testing/testing_protein_domains/sample_orfs_and_hmmscan_df/sample_hmmscan_df.csv'
        sample_hmmscan_df = pd.read_csv(sample_hmmscan_df_path)
        sample_orfs = './test/data_for_testing/testing_protein_domains/sample_orfs_and_hmmscan_df/prodigal_prtn_seq_output.faa'
        converted_df_output = protein_domains.convert_hmmscan_domains_from_protein_to_nuc_pos(sample_hmmscan_df, sample_orfs)

        # Check df
        self.assertEqual(7, converted_df_output.shape[1])
        self.assertEqual(['band', 'contig', 'target_var', 'target_label', 'start', 'stop', 'color'], list(converted_df_output.columns.values))


class test_get_prtn_domain_band_data(unittest.TestCase):

    def test_get_prtn_domain_band_data(self):
        sample_orfs = './test/data_for_testing/testing_protein_domains/sample_orfs_and_hmmscan_df/prodigal_prtn_seq_output.faa'
        hmm_db = '/mnt/nfs/sharknado/Sandbox/Tylo/databases/hmmer/hmmscan/Pfam-A.hmm'
        output_dir = './test/data_for_testing/testing_protein_domains/temp_hmmscan_output'


        actual_df = protein_domains.get_prtn_domain_band_data(sample_orfs, hmm_db, output_dir, 0.01, 'testing_get_prtn_domain_band_data')

        self.assertEqual(7, actual_df.shape[1])
        self.assertEqual(['band', 'contig', 'target_var', 'target_label', 'start', 'stop', 'color'], list(actual_df.columns.values))


        # remove Files
        tblout_path = output_dir + '/testing_get_prtn_domain_band_data_tblout.txt'
        domtblout_path = output_dir + '/testing_get_prtn_domain_band_data_domtblout.txt'
        out_file_path = output_dir + '/testing_get_prtn_domain_band_data_out.txt'
        os.remove(tblout_path)
        os.remove(domtblout_path)
        os.remove(out_file_path)


class test_random_rgb_color(unittest.TestCase):

    def test_random_rgb_color(self):
        result = protein_domains.random_rgb_color()
        self.assertTrue(isinstance(result, str))
        self.assertEqual(result[:4], 'rgb(')


if __name__ == "__main__":
    unittest.main()
