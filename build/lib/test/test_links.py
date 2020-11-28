from fosvis import links

import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import filecmp
import os
import pandas as pd
from pandas.testing import assert_frame_equal
import difflib


class test_get_link_data(unittest.TestCase):

    def test_get_link_data_blastn(self):
        fake_proj_dir = './test/data_for_testing/links_testing/sample_proj_dir'
        list_of_correct_len_contigs_paths = os.listdir('./test/data_for_testing/links_testing/sample_proj_dir/intermediate_outputs/correct_length_contigs')
        list_of_correct_len_contigs_paths = [fake_proj_dir + '/intermediate_outputs/correct_length_contigs/' + s for s in list_of_correct_len_contigs_paths]

        actual_df = links.get_link_data(list_of_correct_len_contigs_paths, 'blastn', fake_proj_dir, 80, 100, '0.60', 50, False)

        # check shape is right
        self.assertEqual(actual_df.shape, (48,7))

        # check column is right
        actual_columns = list(actual_df.columns.values)
        expected_columns = ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos', 'color_thickness']
        self.assertEqual(actual_columns, expected_columns)

    def test_get_link_data_tblastx(self):
        fake_proj_dir = './test/data_for_testing/links_testing/sample_proj_dir'
        list_of_correct_len_contigs_paths = os.listdir('./test/data_for_testing/links_testing/sample_proj_dir/intermediate_outputs/correct_length_contigs')
        list_of_correct_len_contigs_paths = [fake_proj_dir + '/intermediate_outputs/correct_length_contigs/' + s for s in list_of_correct_len_contigs_paths]

        actual_df = links.get_link_data(list_of_correct_len_contigs_paths, 'tblastx', fake_proj_dir, 80, 100, '0.60', 50, False)

        # check shape is right
        self.assertEqual(actual_df.shape, (371,7))

        # check column is right
        actual_columns = list(actual_df.columns.values)
        expected_columns = ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos', 'color_thickness']
        self.assertEqual(actual_columns, expected_columns)


class test_get_pair_links(unittest.TestCase):

    def test_get_pair_blastn_links(self):
        fake_proj_dir = './test/data_for_testing/links_testing/sample_proj_dir'
        query_path = './test/data_for_testing/links_testing/test_get_pair_links/E01.fasta'
        subject_path =  './test/data_for_testing/links_testing/test_get_pair_links/F01.fasta'
        actual_df = links.get_pair_blastn_links(query_path, subject_path, fake_proj_dir, 50, 50, False)

        # check shape is right
        self.assertEqual(actual_df.shape, (16,6))

        # check column is right
        actual_columns = list(actual_df.columns.values)
        expected_columns = ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos']
        self.assertEqual(actual_columns, expected_columns)


    def test_get_pair_tblastx_links(self):
        fake_proj_dir = './test/data_for_testing/links_testing/sample_proj_dir'
        query_path = './test/data_for_testing/links_testing/test_get_pair_links/query.fasta'
        subject_path =  './test/data_for_testing/links_testing/test_get_pair_links/subject.fasta'
        actual_df = links.get_pair_tblastx_links(query_path, subject_path, fake_proj_dir, 50, 50, False)

        # check shape is right
        self.assertEqual(actual_df.shape, (15,6))

        # check column is right
        actual_columns = list(actual_df.columns.values)
        expected_columns = ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos']
        self.assertEqual(actual_columns, expected_columns)


class test_color_similar_links(unittest.TestCase):

    def test_color_similar_links_shape_and_columns(self):

        # load sample data - csv that repersents a df with homology_link_data
        sample_homology_links_df = pd.read_csv('./test/data_for_testing/links_testing/homology_links_pre_color_similar_links_df.csv')

        actual_df = links.color_similar_links(sample_homology_links_df, '0.60', 50)

        # check shape is right
        self.assertEqual(actual_df.shape[1], 7)
        self.assertTrue(actual_df.shape[0] > 0)

        # check column is right
        actual_columns = list(actual_df.columns.values)
        expected_columns = ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos', 'color_thickness']
        self.assertEqual(actual_columns, expected_columns)


class test_merge_lists_with_common_elements(unittest.TestCase):

    def test_merge_lists_with_common_elements(self):
        # Most cases done

        # One item
        listoflist = [[1,2,3]]
        result = links.merge_lists_with_common_elements(listoflist)
        expected_result = [{1,2,3}]
        self.assertEqual(result, expected_result)

        # 4 items - none merged
        listoflist = [[1,2], [3,4], [5,6], [7,8]]
        result = links.merge_lists_with_common_elements(listoflist)
        expected_result = [{1,2}, {3,4}, {5,6}, {7,8}]
        self.assertEqual(result, expected_result)

        # 4 items - all merged
        listoflist = [[1,2,3,4], [4,5,6,7,], [7,8,9,10], [10,11]]
        result = links.merge_lists_with_common_elements(listoflist)
        expected_result = [{1,2,3,4,5,6,7,8,9,10,11}]
        self.assertEqual(result, expected_result)

        # Two items being merged
        listoflist = [[1,2,3], [7,3,9]]
        result = links.merge_lists_with_common_elements(listoflist)
        expected_result = [{1,2,3,7,3,9}]
        self.assertEqual(result, expected_result)

        # Two items being merged, one not
        listoflist = [[1,2,3], [7,3,9], [11,22,33]]
        result = links.merge_lists_with_common_elements(listoflist)
        expected_result = [{1,2,3,7,3,9}, {11,22,33}]
        self.assertEqual(result, expected_result)

        # Multiple merges
        listoflist = [[1,2,3], [7,3,9], [11,22,33], [7,77,777]]
        result = links.merge_lists_with_common_elements(listoflist)
        expected_result = [{1,2,3,7,3,9,7,77,777}, {11,22,33}]
        self.assertEqual(result, expected_result)


class test_calc_percentage_overlap(unittest.TestCase):

    def test_calc_percentage_overlap(self):

        # No overlap
        range1 = (0,10)
        range2 = (10,20)
        result = links.calc_percentage_overlap(range1, range2)
        self.assertEqual(result, 0.0)

        # Normal partial overlap
        range1 = (0,10)
        range2 = (5,15)
        result = links.calc_percentage_overlap(range1, range2)
        self.assertEqual(result, 33.33333333333333)

        # # Normal partial overlap - range1 inverted
        range1 = (10,0)
        range2 = (5,15)
        result = links.calc_percentage_overlap(range1, range2)
        self.assertEqual(result, 33.33333333333333)

        # Normal partial overlap - range1 and range2 inverted
        range1 = (10,0)
        range2 = (15,5)
        result = links.calc_percentage_overlap(range1, range2)
        self.assertEqual(result, 33.33333333333333)

        # Normal partial overlap - range2 inverted
        range1 = (0,10)
        range2 = (15,5)
        result = links.calc_percentage_overlap(range1, range2)
        self.assertEqual(result, 33.33333333333333)

        # Complete overlap
        range1 = (0,100)
        range2 = (10,20)
        result = links.calc_percentage_overlap(range1, range2)
        self.assertEqual(result, 10)


class test_create_color_palette(unittest.TestCase):

    def test_create_color_palette(self):

        # 1 color
        result_1 = links.create_color_palette(1)
        expected_result_1 = ['219.29999999999998,94.65599999999999,86.69999999999999']
        self.assertEqual(result_1, expected_result_1)
        self.assertTrue(all(isinstance(n, str) for n in result_1))

        # 3 color
        result_3 = links.create_color_palette(3)
        expected_result_3 = ['219.29999999999998,94.65599999999999,86.69999999999999',
                             '86.69999999999999,219.29999999999998,94.65599999999999',
                             '94.65599999999999,86.69999999999999,219.29999999999998']
        self.assertEqual(result_3, expected_result_3)
        self.assertTrue(all(isinstance(n, str) for n in result_3))

        # 10 color
        result_10 = links.create_color_palette(10)
        self.assertTrue(all(isinstance(n, str) for n in result_10))





if __name__ == "__main__":
    unittest.main()
