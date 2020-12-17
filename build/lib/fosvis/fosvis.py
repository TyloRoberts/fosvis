from fosvis import seq_util
from fosvis import links
from fosvis import protein_domains
from fosvis import orf
from fosvis import diagram_and_legend

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
from shutil import which

"""
How to use: Script is used through the create_circos_data() and make_diagram() methods.  See documentation pdf for more detail.
"""

############################### create_circos_data #############################

def create_circos_data(contigs, output_dir, project_title, orfs=True, gc=True, custom_histogram=False, custom_histogram_file='', hmm_db='', e_value=0.01, min_contig_length=10000, min_blast_similarity_percentage=90, min_blast_similarity_length=300, link_transperacny='0.60', percent_link_overlap_tolerance=50, include_domains=False, gc_interval_len=100, blast_type='blastn', keep_blast_data=False):
    """
    Create the following input files for a circos diagram:
        - ORF.txt
        - ORF_reverse.txt
        - links.txt
        - karyotype.txt
        - circos.conf

    Args:
        contigs (str): File path to a FASTA file containing contigs
        output_dir (str): File path to a directory where the project folder will be created
        project_title (str): Name of the project directory that is created in the ouput_dir (and used for some file prefixes)
        orfs (bool): If True will include orfs layers, if False don't
        gc (bool): If True include GC layer, if False don't
        custom_histogram (bool): If True, inlcude custom histogram layer, if False don't
        custom_histogram_file (str): File path to a custom histogram data txt file (see documentation for more details)
        hmm_db (str): File path to a pressed hmm database
        e_value (int): The e-value used for the hmmscan search
        min_contig_length (int): The minimum size of a contig to be used in the diagram (idea is that each contig should represent one fosmid)
        min_blast_similarity_percentage (int): Min nucleotide percent similarity (blast percent identity) for link data
        min_blast_similarity_length (int): Min nucleotide similarity length (blast align length) for link data
        link_transperacny (str): The transparency of the links in the diagram in range 0 - 1 (0 is transparent)
        percent_link_overlap_tolerance (int): The amount of overlap percentage required to classify two links as similar enough to color them with the same color
        inlcude_domains (boolean): If True will create protein domain band data, if false will not
        gc_interval_len (int): The interval over which gc content is calculated for the histogram
        blast_type (str): The type of blast to use to create the links (can be 'blastn' OR 'tblastx')
        keep_blast_data (bool): If True will keep the raw blast data, won't if false

    Returns:
        None
    """
    # check if all required programs are in PATH
    if which('circos') is None:
        print("Could not find circos, please install and append to PATH")
        sys.exit()

    if which('prodigal') is None:
        print("Could not find prodigal, please install and append to PATH")
        sys.exit()

    if blast_type == 'blastn':
        if which('blastn') is None:
            print("Could not find blastn, please install and append to PATH")
            sys.exit()

    if blast_type == 'tblastx':
        if which('bl2seq') is None:
            print("Could not find tblastx (called bl2seq in PATH), please install and append to PATH")
            sys.exit()

    if which('hmmscan') is None:
        print("Could not find hmmscan, please install and append to PATH")
        sys.exit()

    project_directory = output_dir + '/' + project_title

    setup_project_dirs(project_directory, blast_type)

    write_paramters_log(project_directory, contigs, output_dir, project_title, hmm_db, e_value, min_contig_length, min_blast_similarity_percentage, min_blast_similarity_length, link_transperacny, percent_link_overlap_tolerance, include_domains, gc, gc_interval_len, blast_type, keep_blast_data)

    # Get correct len contigs
    fosmid_size_contigs = seq_util.remove_too_small_contigs(contigs, min_contig_length)
    correct_size_contigs_path = seq_util.write_correct_len_contigs(project_directory, fosmid_size_contigs)

    # karyotype data
    karyotype_df = seq_util.get_karyotype_data(fosmid_size_contigs)
    karyotype_file = seq_util.write_karyotype_data(project_directory, karyotype_df)

    # Get list of all the correct len contigs
    list_of_correct_len_contigs_paths = os.listdir(project_directory + '/intermediate_outputs/correct_length_contigs')
    list_of_correct_len_contigs_paths.remove('all_correct_size_contigs.fasta')
    list_of_correct_len_contigs_paths = [project_directory + '/intermediate_outputs/correct_length_contigs/' + s for s in list_of_correct_len_contigs_paths]

    # links data
    links_df = links.get_link_data(list_of_correct_len_contigs_paths, blast_type, project_directory, min_blast_similarity_percentage, min_blast_similarity_length, link_transperacny, percent_link_overlap_tolerance, keep_blast_data)
    if links_df.empty:
        print('No links were found, your diagram will still be created but no links will be displayed between your sequences')
    links_df.to_csv(project_directory + '/circos_diagram_input_data/links.txt', sep=' ', index=False, header=False)

    # orf data
    if orfs:
        prodigal_prtn_seq_output = orf.create_orf_data(project_directory, correct_size_contigs_path)

    # Protein domain data
    if include_domains:
        prtn_domain_df = protein_domains.get_prtn_domain_band_data(prodigal_prtn_seq_output, hmm_db, project_directory + '/intermediate_outputs/hmmscan', e_value, project_title)
        prtn_domain_df.to_csv(karyotype_file, mode='a', sep=' ', index=False, header=False)

    # GC data
    if gc:
        gc_content_file = project_directory + '/circos_diagram_input_data/gc_content.txt'
        gc_content_df = seq_util.gc_interval(correct_size_contigs_path, gc_interval_len)
        gc_content_df.to_csv(gc_content_file, sep=' ', index=False, header=False)

    # Create circos conf file
    diagram_and_legend.create_circos_conf_file(project_directory + '/circos_diagram_input_data/circos.conf', project_directory + '/circos_diagram_input_data/karyotype.txt', project_directory + '/circos_diagram_input_data/links.txt', orfs, gc, custom_histogram, orf_forward_file=project_directory + '/circos_diagram_input_data/ORF.txt', orf_reverse_file=project_directory + '/circos_diagram_input_data/ORF_reverse.txt', gc_file=project_directory + '/circos_diagram_input_data/gc_content.txt', custom_histogram_file=custom_histogram_file)
    # if gc:
    #     diagram_and_legend.create_circos_conf_file(project_directory + '/circos_diagram_input_data/circos.conf', project_directory + '/circos_diagram_input_data/karyotype.txt', project_directory + '/circos_diagram_input_data/links.txt', project_directory + '/circos_diagram_input_data/ORF.txt', project_directory + '/circos_diagram_input_data/ORF_reverse.txt', gc_data=project_directory + '/circos_diagram_input_data/gc_content.txt')
    # else:
    #     diagram_and_legend.create_circos_conf_file(project_directory + '/circos_diagram_input_data/circos.conf', project_directory + '/circos_diagram_input_data/karyotype.txt', project_directory + '/circos_diagram_input_data/links.txt', project_directory + '/circos_diagram_input_data/ORF.txt', project_directory + '/circos_diagram_input_data/ORF_reverse.txt')


def setup_project_dirs(project_directory, blast_type):
    """
    Sets up project directory for the project.

    Args:
        project_directory (str): A file path to a directory (the directory for the project)
            that does not already exist yet.
        blast_type (str): The type of blast to use to create the links (can be 'blastn' OR 'tblastx')


    Returns:
        None
    """
    print("Setting up project...")
    try:
        os.mkdir(project_directory)
    except Exception as e:
        print(e)
        print("Was not able to create project directory.")
        sys.exit()

    try:
        os.mkdir(project_directory + '/circos_diagram_input_data')
        os.mkdir(project_directory + '/intermediate_outputs')
        os.mkdir(project_directory + '/intermediate_outputs/prodigal')
        os.mkdir(project_directory + '/intermediate_outputs/correct_length_contigs')
        os.mkdir(project_directory + '/intermediate_outputs/hmmscan')
        if (blast_type == 'blastn'):
            os.mkdir(project_directory + '/intermediate_outputs/blastn')
        elif (blast_type == 'tblastx'):
            os.mkdir(project_directory + '/intermediate_outputs/tblastx')
    except Exception as e:
        print(e)
        print("Was not able to set up internal project directory structure.")
        sys.exit()


def write_paramters_log(project_directory, contigs, output_dir, project_title, hmm_db, e_value, min_contig_length, min_blast_similarity_percentage, min_blast_similarity_length, link_transperacny, percent_link_overlap_tolerance, include_domains, gc, gc_interval_len, blast_type, keep_blast_data):
    """
    Writes the paramters passed into create_circos_data() into the paramters_log.txt

    Args:
        project_directory (str): A file path to the project_directory
        contigs (str): File path to a FASTA file containing contigs
        output_dir (str): File path to a directory where the project folder will be created
        project_title (str): Name of the project directory that is created in the ouput_dir (and used for some file prefixes)
        hmm_db (str): File path to a pressed hmm database
        e_value (int): The e-value used for the hmmscan search
        min_contig_length (int): The minimum size of a contig to be used in the diagram (idea is that each contig should represent one fosmid)
        min_blast_similarity_percentage (int): Min nucleotide percent similarity (blast percent identity) for link data
        min_blast_similarity_length (int): Min nucleotide similarity length (blast align length) for link data
        link_transperacny (str): The transparency of the links in the diagram in range 0 - 1 (0 is transparent)
        percent_link_overlap_tolerance (int): The amount of overlap percentage required to classify two links as similar enough to color them with the same color
        inlcude_domains (boolean): If True will create protein domain band data, if false will not
        gc (boolean): If True the diagram will include a track for GC content, if False it will not
        gc_interval_len (int): The interval over which gc content is calculated for the histogram
        blast_type (str): The type of blast to use to create the links (can be 'blastn' OR 'tblastx')
        keep_blast_data (bool): If True will keep the raw blast data, won't if false

    Returns:
        None
    """
    paramters_log_file_path = project_directory + '/paramaters_log.txt'
    with open(paramters_log_file_path, "x") as f:
        f.write("contigs: " + contigs + '\n'
                "output_dir: " + output_dir + '\n'
                "project_title: " + project_title + '\n'
                "hmm_db: " + hmm_db + '\n'
                "e_value: " + str(e_value) + '\n'
                "min_contig_length: " + str(min_contig_length) + '\n'
                "min_blast_similarity_percentage: " + str(min_blast_similarity_percentage) + '\n'
                "min_blast_similarity_length: " + str(min_blast_similarity_length) + '\n'
                "link_transperacny: " + link_transperacny + '\n'
                "percent_link_overlap_tolerance: " + str(percent_link_overlap_tolerance) + '\n'
                "include_domains: " + str(include_domains) + '\n'
                "gc: " + str(gc) + '\n'
                "gc_interval_len: " + str(gc_interval_len) + '\n'
                "blast_type: " + blast_type + '\n'
                "keep_blast_data: " + str(keep_blast_data) + '\n')
        f.close()


#################### make_diagram ##############################################

def make_diagram(data_dir, ncol=2):
    """
    Makes call to circos to make the diagram.
    Creates the protein domain and karyotype legends.

    Args:
        data_dir (str): File path to a directory containing the following files:
            -> ORF.txt
            -> ORF_reverse.txt
            -> links.txt
            -> karyotype.txt
            -> circos.conf
        ncol (int): Number of columns in the protein_domain legend

    Returns:
        None

    Notes:
        -> Running with -nonparanoid flag so band overlap errors won't terminate the program
    """
    print("Creating circos diagram...")

    circos_stdout_and_stderr_log_path = os.path.dirname(data_dir) + '/intermediate_outputs/circos_stdout_and_stderr_log.txt'
    circos_command = "circos -conf " + data_dir + "/circos.conf -outputdir " + data_dir + " -outputfile circos_diagram.png -noparanoid"

    with open(circos_stdout_and_stderr_log_path,"wb") as stdout:
        subprocess.call([str(circos_command)], stdout=stdout, stderr=subprocess.STDOUT, shell=True)

    diagram_and_legend.prtn_domain_legend(data_dir + "/karyotype.txt", data_dir + "/protein_domain_legend.png", ncol)
    diagram_and_legend.karyotype_legend(data_dir + "/karyotype.txt", data_dir + "/karyotype_legend.png")
