""""
All functions related to link data
"""

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqUtils import GC
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import itertools
import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import statistics
import seaborn as sns
from tqdm import tqdm


def get_link_data(list_of_correct_len_contigs_paths, blast_type, project_directory, min_blast_similarity_percentage, min_blast_similarity_length, link_transperacny, percent_link_overlap_tolerance, keep_blast_data):

    """
    Returns a dataframe containing all the seq similarity link information.
    Colors links that overlap >= percent_link_overlap_tolerance with the same color

    Args:
        list_of_correct_len_contigs_paths list<str>: A list of correct length contig file names (full path, not just file name)
        blast_type (str): The type of blast to be used to create the links (i.e. 'blastn' OR 'tblastx')
        project_directory (str): File path to project directory
        min_blast_similarity_percentage (int): Only inlcude nucleotide similarity for link data with a percent identiity (similarity) >= to this
        min_blast_similarity_length (int): Only include nucleotide similarity for link data with an align_length (length of similarity) >= to this variable
        link_transperacny (str): The transperancy of the links in the diagram in range 0 - 1 (0 is transperant)
        percent_link_overlap_tolerance (int): The amount of overlap percentage required to classify two links as similar enough to color them with the same color
        keep_blast_data (bool): If True will keep the raw blast data, won't if false

    Returns:
        pandas.dataframe: A dataframe with nucleotide homology link data in the form:
            ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos', 'color_thickness']

    Notes:
        - The start position for the origin and terminus does not nessearily have to be bigger than the end pos
    """
    print("Creating Link data...")

    if (blast_type != 'blastn') and (blast_type != 'tblastx'):
        sys.exit("The blast type: " + blast_type  + "is not 'blastn' or 'tblastx'.  It needs to be one of these two.")

    all_homology_links_df = pd.DataFrame(index=[], columns=[])

    # progress bar
    total_len = sum(1 for e in itertools.combinations(list_of_correct_len_contigs_paths, 2))
    with tqdm(total=total_len) as pbar:
        # Loop through all combinations of contigs
        for pair in tqdm(itertools.combinations(list_of_correct_len_contigs_paths, 2)):
            if blast_type == 'blastn':
                pair_homology_df = get_pair_blastn_links(pair[0], pair[1], project_directory, min_blast_similarity_percentage, min_blast_similarity_length, keep_blast_data)
            elif blast_type == 'tblastx':
                pair_homology_df = get_pair_tblastx_links(pair[0], pair[1], project_directory, min_blast_similarity_percentage, min_blast_similarity_length, keep_blast_data)

            # Add new links to main df
            df_to_concat = [all_homology_links_df, pair_homology_df]
            all_homology_links_df = pd.concat(df_to_concat)

            # progress bar
            pbar.update(1)

    # Reset Indexes
    all_homology_links_df = all_homology_links_df.reset_index(drop=True)

    # adjust for circos_indexing
    all_homology_links_df['origin_end_pos'] += 1
    all_homology_links_df['terminus_end_pos'] += 1

    return color_similar_links(all_homology_links_df, link_transperacny, percent_link_overlap_tolerance)


def get_pair_blastn_links(query, subject, project_directory, min_blast_similarity_percentage, min_blast_similarity_length, keep_blast_data):
    """
    Indentify regions of nucleotide homology between two SeqRecord objects (greater
    than 'min_blast_similarity_percentage' across intervals of more than
    'min_blast_similarity_length' bp) and return a dataframe containing the infomration.

    Args:
        query (str): File path to a contig fasta file to be the query
        subject (str): File path to a contig fasta file to be the subject
        project_directory (str): File path for the project directory
        min_blast_similarity_percentage (int): Only use results with a percent identiity >= to this variable
        min_blast_similarity_length (int): Only use results that have an align_length >= to this variable
        keep_blast_data (bool): If True will keep the raw blast data, won't if false

    Returns:
        A dataframe conatining the homolgy information with the following structure is returned with each row repersenting a differnt alignment:
            ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos']

    Notes:
    - In the output the returned origin/terminus start positions do not have to be > than the end pos
    """

    # Run BLAST
    blastn_out = project_directory + "/intermediate_outputs/blastn/blastn_" + os.path.basename(query).replace('.fasta', '') + "_to_" + os.path.basename(subject).replace('.fasta', '') + "_output.xml"
    output = NcbiblastnCommandline(query=query, subject=subject, outfmt=5, out=blastn_out)()[0]

    # Read blastn output
    with open(blastn_out, "r") as blast_results_file:
        blast_result_record = NCBIXML.read(blast_results_file)

    seq_homology_df = pd.DataFrame(index=[], columns=['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos'])

    # Parse results
    for alignment in blast_result_record.alignments:
        subject_id = alignment.hit_def.split()[0]
        query_id = blast_result_record.query.split()[0]

        for hsp in alignment.hsps:
            identities_percentage = (hsp.identities / hsp.align_length) * 100

            # Append to df if similar/long enough
            if (identities_percentage >= min_blast_similarity_percentage) and (hsp.align_length >= min_blast_similarity_length):
                new_alignment = {'origin_var' : query_id, 'origin_start_pos' : hsp.query_start, 'origin_end_pos' : hsp.query_end, 'terminus_var' : subject_id, 'terminus_start_pos' : hsp.sbjct_start, 'terminus_end_pos' : hsp.sbjct_end}
                seq_homology_df = seq_homology_df.append(new_alignment, ignore_index=True)

    if keep_blast_data == False:
        os.remove(blastn_out)

    return seq_homology_df


def get_pair_tblastx_links(query, subject, project_directory, min_blast_similarity_percentage, min_blast_similarity_length, keep_blast_data):
    """
    Indentify regions of protein similarity using tblastx between two SeqRecord objects (greater
    than 'min_blast_similarity_percentage' similarity across intervals of more than
    'min_blast_similarity_length' bp) and return a dataframe containing the infomration.

    Args:
        query (str): File path to a contig fasta file to be the query
        subject (str): File path to a contig fasta file to be the subject
        project_directory (str): File path for the project directory
        min_blast_similarity_percentage (int): Only use results with a percent identiity >= to this variable
        min_blast_similarity_length (int): Only use results that have an aligned length >= to this variable
        keep_blast_data (bool): If True will keep the raw blast data, won't if false

    Returns:
        A dataframe conatining the homolgy information with the following structure is returned with each row repersenting a differnt alignment:
            ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos']

    Notes:
    - CHECK THIS In the output the returned origin/terminus start positions do not have to be > than the end pos
    - Alignment length is calculated abs(Query_End - Query_Start)
    """

    # Run Create output
    tblastx_out = project_directory + "/intermediate_outputs/tblastx/tblastx_" + os.path.basename(query) + "_to_" + os.path.basename(subject) + "_output.txt"

    # Run tblastx
    tblastx_command = 'bl2seq -p tblastx -i ' + query + ' -j ' + subject + ' -o ' + tblastx_out
    subprocess.call([str(tblastx_command)], shell=True)

    # Parse output for all matches
    results = parse_tblastx(tblastx_out)

    similarity_filt = results['Identities_Percent'] >= min_blast_similarity_percentage
    length_filt = abs(results['Query_End'] - results['Query_Start']) >= min_blast_similarity_length

    results = results[similarity_filt & length_filt]

    # Remove unnesary columns from results
    results.drop('Score', axis=1, inplace=True)
    results.drop('Expect', axis=1, inplace=True)
    results.drop('Identities_Percent', axis=1, inplace=True)
    results.drop('Positives_Percent', axis=1, inplace=True)

    # Arrange columns in correct order
    results = results[['Query', 'Query_Start', 'Query_End', 'Subject', 'Sbjct_Start', 'Sbjct_End']]

    # Change name of the columns
    results.columns = ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos']

    if keep_blast_data == False:
        os.remove(tblastx_out)

    return results


def parse_tblastx(input):
    """
    Parse a tblastx output into a pandas df output.
    Meant for an ouptut from 2 seqeunces

    Args:
        input (str):
        output (str):

    Returns:
        pandas df: Returns a pandas df with the columns [Query, Subject, Score, Expect, Identities_Percent, Positives_Percent,
                                                         Query_Start, Query_End, Sbjct_Start, Sbjct_End]

    Notes:
    - Score is in bits
    - Some e values don't have a number at beggining (e.g. e-163), a 1 is added to front (e.g. 1e-163)
    """

    with open(input, 'r') as f:
        content = f.readlines()

    # Remove white space
    content = [x.strip() for x in content]

    df_data = {
      'Query': [],
      'Subject': [],
      'Score': [],
      'Expect': [],
      'Identities_Percent': [],
      'Positives_Percent': [],
      'Query_Start': [],
      'Query_End': [],
      'Sbjct_Start': [],
      'Sbjct_End': []
    }

    # Dataframe Row
    Query = str(content[0].replace('Query= ', ''))
    Subject = 0
    Score = 0
    Expect = 0
    Identities_Percent = 0
    Positives_Percent = 0
    Query_Start = 0
    Query_End = 0
    Sbjct_Start = 0
    Sbjct_End = 0

    for line in content[1:]:

        # Get Sbjct seq name
        if line.startswith('>'):
            Subject = str(line.replace('>', ''))

        if line.startswith('Score'):
            # Write row to dataframe - signals at end of a data entry
            # Unless its the first run through
            if Score != '':
                df_data["Query"].append(Query)
                df_data["Subject"].append(Subject)
                df_data["Score"].append(Score)
                df_data["Expect"].append(Expect)
                df_data["Identities_Percent"].append(Identities_Percent)
                df_data["Positives_Percent"].append(Positives_Percent)
                df_data["Query_Start"].append(Query_Start)
                df_data["Query_End"].append(Query_End)
                df_data["Sbjct_Start"].append(Sbjct_Start)
                df_data["Sbjct_End"].append(Sbjct_End)

            # Get Score/Expect Data
            Score = float(line.replace('Score = ', '').split()[0])
            temp_expect = line.split(', ')[1].split()[2]
            if temp_expect.startswith('e'):
                Expect = float('1' + temp_expect)
            else:
                Expect = float(temp_expect)

        if line.startswith('Identities'):
            Identities_Percent = float(line.split('(')[1].split('%')[0])
            Positives_Percent = float(line.split(',')[1].split('(')[1].split('%')[0])

        if line.startswith('Query'):
            # If Query/Sbjct are multiple lines value will just be updated until
            # they are written to the df
            Query_Start = int(line.split()[1])
            Query_End = int(line.split()[3])

        if line.startswith('Sbjct'):
            # If Query/Sbjct are multiple lines value will just be updated until
            # they are written to the df
            Sbjct_Start = int(line.split()[1])
            Sbjct_End = int(line.split()[3])

    df = pd.DataFrame(df_data)

    return df


def color_similar_links(homolgy_links_df, link_transperacny, percent_link_overlap_tolerance):
    """
    Given a set of homology links, finds links that repersent the same seqeunce (they overlap with each other) being repersented by
    differnt links and color those with a matching unique color.

    Args:
        homolgy_links_df (pandas.df): A pandas dataframe repersenting all the nucleotide homology link data for a circos diagram
            ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos']
        link_transperacny (str): The transperancy of the links in the diagram in range 0 - 1 (0 is transperant)
        percent_link_overlap_tolerance (int): The amount of overlap percentage required to classify two links as similar enough to color them with the same color

    Returns:
        (pandas.df): A dataframe with link data colored accoridng to links that repersent the same sequence.
            ['origin_var', 'origin_start_pos', 'origin_end_pos', 'terminus_var', 'terminus_start_pos', 'terminus_end_pos', 'color_thickness']

    Notes:
        - Required that the dataframe indexes are from 0 - n in increasing order (groups are made based on indexes so they must all be distinct)
    """

    # List of link pairs that overlap sufficiently
    matching_link_indexes = []

    # Go thorugh all combinations of links - if overlap sufficiently, append to matching_link_indexes
    for index, link in homolgy_links_df.iterrows():

        for pot_match_index, pot_match_link in homolgy_links_df.iterrows():
            if pot_match_index != index: # don't compare to itself

                origin_percent_overlap = calc_percentage_overlap((link['origin_start_pos'], link['origin_end_pos']), (pot_match_link['origin_start_pos'], pot_match_link['origin_end_pos']))
                terminus_percent_overlap = calc_percentage_overlap((link['origin_start_pos'], link['origin_end_pos']), (pot_match_link['terminus_start_pos'], pot_match_link['terminus_end_pos']))

                # Case: origin of link to origin of pot_match_link
                if (link['origin_var'] == pot_match_link['origin_var']) and (origin_percent_overlap >= percent_link_overlap_tolerance): # same contig and overlap
                    matching_link_indexes.append([index, pot_match_index])

                # Case: origin of link to terminus of pot_match_link
                if (link['origin_var'] == pot_match_link['terminus_var']) and (terminus_percent_overlap >= percent_link_overlap_tolerance): # same contig and overlap
                    matching_link_indexes.append([index, pot_match_index])

        for pot_match_index, pot_match_link in homolgy_links_df.iterrows():
            if pot_match_index != index:

                origin_percent_overlap = calc_percentage_overlap((link['terminus_start_pos'], link['terminus_end_pos']), (pot_match_link['origin_start_pos'], pot_match_link['origin_end_pos']))
                terminus_percent_overlap = calc_percentage_overlap((link['terminus_start_pos'], link['terminus_end_pos']), (pot_match_link['terminus_start_pos'], pot_match_link['terminus_end_pos']))

                # Case: terminus of link to origin of pot_match_link
                if (link['terminus_var'] == pot_match_link['origin_var']) and (origin_percent_overlap >= percent_link_overlap_tolerance):
                    matching_link_indexes.append([index, pot_match_index])

                # Case: terminus of link to terminus of pot_match_link
                if (link['terminus_var'] == pot_match_link['terminus_var']) and (terminus_percent_overlap >= percent_link_overlap_tolerance):
                    matching_link_indexes.append([index, pot_match_index])


    # Put links that repersent the same sequence in sets
    link_indexs_w_common_homology = merge_lists_with_common_elements(matching_link_indexes)

    # Loop through all indexes in homolgy_links_df, if they are not in link_indexs_w_common_homology
    # then append a set with just that index i.e. {index} to another list and then combine with link_indexs_w_common_homology
    not_in_link_indexs_w_common_homology = []
    for index in list(homolgy_links_df.index.values):

        # Check if index is not in any of the sets in link_indexs_w_common_homology
        if not any(index in sublist for sublist in link_indexs_w_common_homology):
            not_in_link_indexs_w_common_homology.append({index})

    # Contains every index in a set by itself or a set with others
    link_indexs_w_common_homology = link_indexs_w_common_homology + not_in_link_indexs_w_common_homology

    # create seaborn color palette with len of link_indexs_w_common_homology
    color_palette = create_color_palette(len(link_indexs_w_common_homology))

    # Create dict with indexes as keys and colors as values - sets every index to same color as its group (if its in a group)
    link_index_color_dict = {}
    color_palette_counter = 0
    for index_set in link_indexs_w_common_homology:
        color = '(' + color_palette[color_palette_counter] + "," + link_transperacny + ')'
        for link_index in index_set:
            link_index_color_dict[link_index] = color
        color_palette_counter += 1

    # Go thorugh dict and set links to their correct color in homolgy_links_df
    for key, value in link_index_color_dict.items():
        homolgy_links_df.at[key, 'color_thickness'] = 'color=' + value + ',thickness=10p'

    return homolgy_links_df


def merge_lists_with_common_elements(listoflist):
    """
    Given a list of lists, merge all lists that share a common element, and repeat until there are no more lists with the same item.

    Args:
        listoflist (list<list<item>>): A list of lists to be merged

    Returns:
        (list<set>): A list of sets of merged items

    Source:
        Based off of https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements
    """


    out = []
    while len(listoflist)>0:
        first, *rest = listoflist
        first = set(first)

        lf = -1
        while len(first)>lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2

        out.append(first)
        listoflist = rest

    return out


def calc_percentage_overlap(range1, range2):
    """
    Calculates the percentage of overlap betwen two ranges.

    Args:
        range1 (tuple<int>): A range repersented as a tuple
        range2 (tuple<int>): A second range repersented as a tuple

    Returns:
        (int): The percentage overlap between ranges A and B

    Notes:
    - Ranges don't have to be passed in a smalles to largest value
    - The percentage of overlap between the two ranges is calculated by the
      following (simplified) formula:
        percentage_overlap = (width of overlaping region / total range of both ranges combined) * 100
    """

    A = range1
    B = range2

    if range1[0] > range1[1]:
        A = (range1[1], range1[0])

    if range2[0] > range2[1]:
        B = (range2[1], range2[0])

    return abs((statistics.median([A[0],B[0],B[1]])-statistics.median([A[1],B[0],B[1]])) / (max(A[1],B[1])-min(A[0],B[0]))) * 100


def create_color_palette(num_colors):
    """
    Creates a color palette with a specified amount of rgb colors.
    Uses the seaborn hls color palette

    Args:
        num_colors (int): The number of colors to be returned

    Returns:
        list<str>: A list of rgb tuples in the form "<num1>,<num2>,<num3>"
    """
    f = lambda y : tuple(map(lambda x : x*255, y))
    color_palette = [f(color) for color in sns.color_palette("hls", num_colors)]

    str_color_palette = []
    for color in color_palette:
        str_color_palette.append(",".join(map(str,color)))

    return str_color_palette
