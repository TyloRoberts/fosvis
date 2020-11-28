"""
Functions that work with seqeunces.
Realtes to karyotype and GC data
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

def remove_too_small_contigs(contigs, min_contig_length):
    """
    Loop thorugh a fasta file and return a list of BioPython SeqRecord objects that repersent all the seqeunces
    in the fasta file with length >= min_contig_length

    Args:
        contigs (str): A file path to a .fasta file containing contigs
        min_contig_length (int): Contigs smaller than this size won't be added to the returned list

    Returns:
        A list of BioPython SeqRecord objects that repersent only the seqeunces in the contigs file with
        size >= min_contig_length
    """

    print("Removing too small contigs...")
    num_contigs_orginally = 0
    correct_size_contigs = []
    for record in SeqIO.parse(contigs, "fasta"):
        if len(record) >= min_contig_length:
            correct_size_contigs.append(record)
        num_contigs_orginally += 1

    if len(correct_size_contigs) == 0:
        sys.exit("No contigs with length >= min_contig_length were found.")

    print(str(len(correct_size_contigs)) + '/' + str(num_contigs_orginally) + ' contigs were of correct size')

    return correct_size_contigs


def write_correct_len_contigs(project_directory, fosmid_size_contigs):
    """
    Writes all the fosmid_size_contigs to one file and returns that file path and
    writes each seqeunce in fosmid_size_contigs to a file individually

    Args:
        project_directory (str): A file path to the project_directory
        fosmid_size_contigs (list<Biopython.SeqRecord>): A list of BioPython SeqRecord objects that repersent only the seqeunces in the contigs file with
            size >= min_contig_length

    Return:
        correct_size_contigs_path (str): A file path to the file containing all the fosmid_size_contigs together

    """
    # Write all correct size contigs to their own fasta file
    for seq in fosmid_size_contigs:
        seq_path = project_directory + "/intermediate_outputs/correct_length_contigs/" + seq.id + ".fasta"
        SeqIO.write(seq, seq_path, "fasta")

    # writes all the correct_size_contigs to one file
    correct_size_contigs_path = project_directory + "/intermediate_outputs/correct_length_contigs/all_correct_size_contigs.fasta"
    all_seqs = []
    for seq in fosmid_size_contigs:
        all_seqs.append(seq)
    SeqIO.write(all_seqs, correct_size_contigs_path, "fasta")

    return correct_size_contigs_path


def get_karyotype_data(fosmid_size_contigs):
    """
    Creates the karyotype data for the sequences in fosmid_size_contigs in the form of a dataframe.

    Args:
        fosmid_size_contigs (list<SeqRecord>): A list of BioPython SeqRecord objects that the karyotype is to be made for

    Returns:
        karyotype_df (pandas.dataframe): A dataframe containing the karyotype data in the df:
            ['chr_prefix', '-prefix', 'variable_name', 'diagram_label', 'start', 'end', 'color']

    Notes:
    -> chr_prefix col is just 'chr'
    -> -prefix col is just '-'
    -> variable_name is the contig
    -> Diagram label is a number from 1 - # of contigs
    -> start is 1
    -> end is (length + 1) because of circos indexing (i.e. bp 1 is repersented by the range 1 - 2)
    """

    print("Creating karyotype data...")
    karyotype_df = pd.DataFrame(index=[], columns=['chr_prefix', '-prefix', 'variable_name', 'diagram_label', 'start', 'end', 'color'])

    counter = 1
    for seq in fosmid_size_contigs:
        karyotype_row = {'chr_prefix' : 'chr', '-prefix' : '-', 'variable_name' : seq.id, 'diagram_label' : str(counter), 'start' : 1, 'end' : len(seq) + 1, 'color' : 'rgb(120,120,120,0.4)'}
        karyotype_df = karyotype_df.append(karyotype_row, ignore_index=True)
        counter += 1

    return karyotype_df


def write_karyotype_data(project_directory, karyotype_df):
    """
    Writes a df containing karyotype data to a txt file in the project_directory

    Args:
        project_directory (str): A file path to the project_directory
        karyotype_df (pandas.DataFrame): A dataframe contaiing karyotype data

    Returns:
        karyotype_file (str): The file path to where the df was written

    """
    karyotype_file = project_directory + '/circos_diagram_input_data/karyotype.txt'
    karyotype_df.to_csv(karyotype_file, sep=' ', index=False, header=False)
    return karyotype_file


def gc_interval(seqs_file, interval_length):
    """
    Returns a dataframe with the GC content of a seqeunce in interval lengths

    Args:
        seqs_file (str): A file with contigs to determine GC content for (can contain multiple seqs)
        interval_length (int): The length of the interval to determine the GC-content in

    Returns:
        (pandas.core.frame.DataFrame): A dataframe with GC-content over set intervals of a seqeunce in the form:
            ['contig', 'interval_start', 'interval_end', 'gc_content']

    Notes:
    - Bipython seq[start:end] is start (0-based) to end (0-based, non-inlcusive)
    - Due to circos range indexing the GC interval for bp pos x - y would be x - (y + 1)
        - The + 1 is so it covers full range of last bp
        - e.g. gc content for bp's 5 - 10 would be ['contig', 5, 11, 'gc_content']
    """
    print("Creating GC content data...")
    gc_content_df = pd.DataFrame(columns = ['contig', 'interval_start', 'interval_end', 'gc_content'], index = [])

    num_seqs = len([1 for line in open(seqs_file) if line.startswith(">")])
    with tqdm(total=num_seqs) as pbar:
        for sequence in SeqIO.parse(seqs_file, 'fasta'):

            current_interval_start = 0

            while ((current_interval_start + interval_length) < (len(sequence.seq))):
                seq_subset = sequence.seq[current_interval_start:current_interval_start+interval_length]
                gc_content = GC(seq_subset)
                gc_new_row = {'contig' : sequence.id,
                              'interval_start': current_interval_start + 1,
                              'interval_end': (current_interval_start + (interval_length + 1)), 'gc_content': gc_content} # good
                gc_content_df = gc_content_df.append(gc_new_row, ignore_index=True)

                current_interval_start += (interval_length)

            else:
                seq_subset = sequence.seq[current_interval_start:]
                gc_content = GC(seq_subset)
                gc_new_row = {'contig' : sequence.id,
                              'interval_start': current_interval_start + 1,
                              'interval_end': len(sequence.seq) + 1, 'gc_content': gc_content}
                gc_content_df = gc_content_df.append(gc_new_row, ignore_index=True)

            pbar.update(1)


    return gc_content_df
