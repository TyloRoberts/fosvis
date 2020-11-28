"""All functions related to ORF data"""

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



def create_orf_data(project_directory, correct_size_contigs_path):
    """
    Calls the functions to run prodigal, parse the output and write the results to
    txt files.

    Args:
        project_directory (str): A file path to the project_directory
        correct_size_contigs_path (str): A file path to the file containing all the correct sized contigs together

    Returns:
        prodigal_prtn_seq_output (str): The file path to prodigal's prtn translation output (.faa file)
    """

    print("Creating ORF data...")

    prodigal_output = project_directory + '/intermediate_outputs/prodigal/prodigal_orf_output'
    prodigal_prtn_seq_output = project_directory + '/intermediate_outputs/prodigal/prodigal_prtn_seq_output.faa'
    run_prodigal(correct_size_contigs_path, prodigal_output, prodigal_prtn_seq_output)

    prodigal_forward_and_reverse_df = parse_prodigal_sco_output(prodigal_output)

    prodigal_forward_and_reverse_df[0].to_csv(project_directory + '/circos_diagram_input_data/ORF.txt', sep=' ', index=False, header=False)
    prodigal_forward_and_reverse_df[1].to_csv(project_directory + '/circos_diagram_input_data/ORF_reverse.txt', sep=' ', index=False, header=False)

    return prodigal_prtn_seq_output


def run_prodigal(input, output, prtn_output, format='sco'):
    """
    Runs prodigal on the input .fasta file

    Args:
        input (str): File path to an input file (.fasta) to run prodigal on
        output (str): File path to output file (in format of format paramter)
        prtn_output (str): File path to a output file contianing protein translations (should be a .faa file)
        format (str): output format of the output file (default = sco)
            gbk:  Genbank-like format
            gff:  GFF format
            sqn:  Sequin feature table format
            sco:  Simple coordinate output

    Returns:
        None

    Notes:
        -c flag (Closed ends.  Do not allow genes to run off edges.) was not used.
        - When the ORF is indicated as -1 (or reverse) this means the orf was found on the reverse complement,
          in this case the indicated start is actually where there gene would have ended (just going in the other direction) and the
          end is where it would have started relative to the input sequence
        - Prodigal uses GTG as M (start codon)
        - Sequence must be 20000 characters for prodigal to run
    """
    prodigal_stdout_log_path = os.path.dirname(output) + '/prodigal_stdout_and_stderr_log.txt'

    prodigal_command = "prodigal -i " + input + " -o " + output + " -a " + prtn_output + " -f " + format
    with open(prodigal_stdout_log_path,"wb") as stdout:
        subprocess.call([str(prodigal_command)], stdout=stdout, stderr=subprocess.STDOUT, shell=True)

def parse_prodigal_sco_output(prodigal_sco_output):
    """
    Parses prokka sco format output file into a pandas DataFrame

    Args:
        prodigal_sco_output (str): File path to a prokka sco format ouptut file

    Returns:
        Two pandas DataFrames (forward strand genes, reverse strand genes)
            ['contig', 'orf_start', 'orf_end']

    Notes:
    - For reverse_strand ORFs, the orf_start is really the end (because it is in reverse)
    """

    with open(prodigal_sco_output, 'r') as f:

        forward_strand_orf_df = pd.DataFrame(index=[], columns=['contig', 'orf_start', 'orf_end'])
        reverse_strand_orf_df = pd.DataFrame(index=[], columns=['contig', 'orf_start', 'orf_end'])
        curContig = ""

        for line in f:
            if (line.startswith("# Sequence Data")):
                curContig = line.split('seqhdr="', 1)[1].rstrip()[:-1].split()[0]

            if (line.startswith(">")):
                start = line.split("_")[1]
                end = line.split("_")[2]
                strand = line.split("_")[3].rstrip()

                if strand == '+':
                    # add 1 to orf_end to adjust for circos indexing
                    new_forward_orf = {'contig' : curContig, 'orf_start' : int(start), 'orf_end' : int(end) + 1}
                    forward_strand_orf_df = forward_strand_orf_df.append(new_forward_orf, ignore_index=True)

                if strand == '-':
                    # add 1 to orf_end (really the start of the reverse strand orf) to adjust for circos indexing
                    new_reverse_orf = {'contig' : curContig, 'orf_start' : int(start), 'orf_end' : int(end) + 1}
                    reverse_strand_orf_df = reverse_strand_orf_df.append(new_reverse_orf, ignore_index=True)

        # get all values as ints
        forward_strand_orf_df['orf_start'] = forward_strand_orf_df['orf_start'].astype(int)
        forward_strand_orf_df['orf_end'] = forward_strand_orf_df['orf_end'].astype(int)
        reverse_strand_orf_df['orf_start'] = reverse_strand_orf_df['orf_start'].astype(int)
        reverse_strand_orf_df['orf_end'] = reverse_strand_orf_df['orf_end'].astype(int)

        return forward_strand_orf_df, reverse_strand_orf_df
