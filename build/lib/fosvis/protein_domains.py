"""
All functions related to protein_domain data and HMMER
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


def run_hmmscan(orfs, hmm_db, output_dir, E_value, prefix):
    """
    Runs hmmscan on open reading frames (ORFs) using a specified (pressed) HMM database

    Args:
        orfs (str): File path to .faa file containing all the ORFs
        hmm_db (str): File path to a pressed HMM database
        output_dir (str): Directory to to put all output files in (should already exist)
        E_value (int): The E value used when running hmmscan
        prefix (str): Prefix for all the output files (e.g. the sample name)

    Returns:
        None
    """

    out = output_dir + '/' + prefix + '_out.txt'
    tblout = output_dir + '/' + prefix + '_tblout.txt'
    domtblout = output_dir + '/' + prefix + '_domtblout.txt'
    hmmscan_stdout_and_stderr_log_path = output_dir + '/hmmscan_stdout_and_stderr_log.txt'

    hmmscan_command = 'hmmscan -o ' + out + ' --tblout ' + tblout + ' --domtblout ' + domtblout + ' -E ' + str(E_value) + ' ' + hmm_db + ' ' + orfs

    with open(hmmscan_stdout_and_stderr_log_path, "wb") as stdout:
        subprocess.call([str(hmmscan_command)], stdout=stdout, stderr=subprocess.STDOUT, shell=True)


def parse_hmmscan_domtblout(domtblout):
    """
    Parse a hmmscan domtbl output file and return a pandas dataframe with the output.
    Designed to be used for a hmmscan output file that was for hmmscan run on a prodigal output.

    Args:
        domtblout (str): File path to a hmmscan domtbl output file to be parsed

    Returns:
        (pandas df): Contains the parsed domtblout information, in the form:
            ['target', 'orf_id', 'orf_start', 'orf_stop', 'target_description']
    """

    with open(domtblout, 'r') as f:

        hmmscan_domains = pd.DataFrame(index=[], columns=['target', 'orf_id', 'orf_start', 'orf_stop', 'target_description'])

        for line in f:

            if (line[0] == '#'):
                continue

            target = line.split()[0]
            orf_id = line.split()[3]
            orf_start = line.split()[17]
            orf_stop = line.split()[18]
            target_description = ' '.join(line.split()[22:])

            new_target_entry = {'target' : target, 'orf_id': orf_id, 'orf_start' : orf_start, 'orf_stop' : orf_stop, 'target_description' : target_description}
            hmmscan_domains = hmmscan_domains.append(new_target_entry, ignore_index=True)

    return hmmscan_domains


def convert_hmmscan_domains_from_protein_to_nuc_pos(hmmscan_df, orf_file):
    """
    Takes in a hmmscan dataframe and determines the nucleotide start/stop positions
    of the idneitfied protein domain in terms of its protein start/stop pos given by the
    hmmscan output.

    Also colors the protein domains the same if they are the same domain.

    Args:
        hmmscan_df (pandas.dataframe): A dataframe from a parsed hmmscan domtbl file
            Needs to be in the following format: ['target', 'orf_id', 'orf_start', 'orf_stop', 'target_description']
        orf_file (str): File path to prodigal orf protein seqeunce file
            Headers need to be in the following format (format from prodigal):
            <contig>_<orf number identifed> # <bp start> # <bp stop> # <strand> # etc.

    Returns:
        (pandas.dataframe): A dataframe in the following fomrat ['band', 'contig', 'target_var', 'target_label', 'start', 'stop', 'color']
    """

    # Create dictionary of orf_id and [start, stop, strand]
    orf_start_stop_dict = {}
    for orf in SeqIO.parse(orf_file, 'fasta'):
        orf_start = orf.description.split()[2]
        orf_stop = orf.description.split()[4]
        forward_or_reverse = orf.description.split()[6] # 1 or -1
        orf_start_stop_dict [orf.id] = [orf_start, orf_stop, forward_or_reverse]

    band_data_df = pd.DataFrame(index=[], columns=['band', 'contig', 'target_var', 'target_label', 'start', 'stop', 'color'])

    # Make Dictionary with keys of every unique domain and values as a random rgb color
    # loop thorugh every row of dataframe and make list of all the domains and turn that into a set
    unique_domains = {}
    for index, row in hmmscan_df.iterrows():
        random_color = random_rgb_color()
        unique_domains[row['target']] = random_color


    for index, row in hmmscan_df.iterrows():

        # Get the unique color for that specific target (all bands with that target will be same color)
        band_color = unique_domains[row['target']]

        contig = row['orf_id'].rsplit('_', 1)[0]

        orf_bp_start = orf_start_stop_dict[row['orf_id']][0]
        orf_bp_stop = orf_start_stop_dict[row['orf_id']][1]
        orf_forward_or_reverse = orf_start_stop_dict[row['orf_id']][2]

        domain_bp_start = None
        domain_bp_stop = None

        # ORF on forward strand
        if int(orf_forward_or_reverse) == 1:
            domain_bp_start = int(orf_bp_start) + 3 * (int(row['orf_start']) - 1)
            domain_bp_stop = int(orf_bp_start) + 3 * int(row['orf_stop']) - 1
        # ORF on reverse strand
        elif int(orf_forward_or_reverse) == -1:
            domain_bp_stop = int(orf_bp_stop) - 3 * (int(row['orf_start']) - 1)
            domain_bp_start = int(orf_bp_stop) - 3 * int(row['orf_stop']) + 1

        target_description_spaces_to_underscore = row['target_description'].replace(' ', '_')

        new_band = {'band' : 'band', 'contig' : contig, 'target_var' : row['target'], 'target_label' : target_description_spaces_to_underscore, 'start' : str(domain_bp_start), 'stop' : str(domain_bp_stop), 'color' : band_color}
        band_data_df = band_data_df.append(new_band, ignore_index=True)

    return band_data_df


def get_prtn_domain_band_data(orfs, hmm_db, output_dir, E_value, prefix):
    """
    Calls functions to run hmmscan, parse the output, convert the start/stop positions
    of the domains from positions on the ORFs to positions on the original contig and
    creates a dataframe with the data in the correct format for circos band data.
    This function is meant to be used with ORFs with headers created by the prodigal tool.

    Args:
        orfs (str): File path to .faa file containing all the ORFs (should be created by prodigal and in the prodigal format)
        hmm_db (str): File path to a pressed HMM database
        output_dir (str): Directory to to put all hmmscan output files in (should already exist)
        E_value (str): The E value used when running hmmscan
        prefix (str): Prefix for all the output files (e.g. the sample name)

    Returns:
        (pandas.dataframe): A dataframe with the band data in the form:
            ['band', 'contig', 'target_var', 'target_label', 'start', 'stop', 'color']
    """
    print("Creating protein domain band data - this may take a while depending on amount/size of fosmids...")
    run_hmmscan(orfs, hmm_db, output_dir, E_value, prefix)
    parsed_hmm_df = parse_hmmscan_domtblout(output_dir + "/" + prefix + "_domtblout.txt")
    return convert_hmmscan_domains_from_protein_to_nuc_pos(parsed_hmm_df, orfs)


def random_rgb_color():
    """
    Generate a random RGB Color.

    Returns:
        (str): A string in the fomrat 'rgb(x,y,z)'
    """

    color = list(np.random.choice(range(256), size=3))
    return 'rgb(' + str(color[0]) + "," + str(color[1]) + ',' + str(color[2]) + ')'
