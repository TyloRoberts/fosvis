"""All functions related to creating diagrams and legends"""

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



def find_layer_radii(orfs, gc, custom_histogram):
    """
    Creates a dictionary with the proper circos radii for the given set of layer inputs to be used in the conf file

    Args:
        orfs (bool): If True will include orfs layers, if False don't
        gc (bool): If True include GC layer, if False don't
        custom_histogram (bool): If True, inlcude custom histogram layer, if False don't

    Returns:
        (tuple(int, dict)): A tuple with the first value being an int that corresponds to a layer combination
                            (as seen in the notes section below) and the second value is a dictioary that corresponds
                            to the correct circos radius positions for that layer combiantion

    Notes:
        Options:
        0 = karyotype, links
        1 = karyotype, ORFs, links
        2 = karyotype, GC, links
        3 = karyotype, custom_histogram, links
        4 = karyotype, ORFs, GC, links
        5 = karyotype, ORFs, custom_histogram, links
        6 = karyotype, GC, custom_histogram, links
        7 = karyotype, ORFs, GC, custom_histogram, links
    """

    # possible combiantions of layer radii
    layer_radii_combinations = [{'karyotype' : '0.80r', 'links' : '0.96r'},
                                {"karyotype": '0.80r', "ORF": ('0.91r', '0.96r'), "ORF_reverse" : ('0.83r', '0.88r'), "links" : '0.79r'},
                                {"karyotype": '0.80r', "GC": ('0.84r', '0.96r'), "links" : '0.80r'},
                                {"karyotype": '0.80r', "custom_histogram" : ('0.84r', '0.96r'), "links" : '0.80r'},
                                {"karyotype": '0.80r', "ORF": ('0.91r', '0.96r'), "ORF_reverse" : ('0.83r', '0.88r'), "GC" :  ('0.68r', '0.8r'), "links" : '0.64r'},
                                {"karyotype": '0.80r', "ORF": ('0.91r', '0.96r'), "ORF_reverse" : ('0.83r', '0.88r'), "custom_histogram" : ('0.68r', '0.8r'), "links" : '0.64r'},
                                {"karyotype": '0.80r', "GC" : ('0.84r', '0.96r'), "custom_histogram" : ('0.68r', '0.80r'), "links" : '0.64r'},
                                {"karyotype": '0.80r', "ORF": ('0.91r', '0.96r'), "ORF_reverse" : ('0.83r', '0.88r'), "GC" : ('0.68r', '0.8r'), "custom_histogram" : ('0.52r', '0.64r'), "links" : '0.48r'}]

    if (orfs == False and gc == False and custom_histogram == False):
        return (0, layer_radii_combinations[0])

    elif (orfs == True and gc == False and custom_histogram == False):
        return (1, layer_radii_combinations[1])

    elif (orfs == False and gc == True and custom_histogram == False):
        return (2, layer_radii_combinations[2])

    elif (orfs == False and gc == False and custom_histogram == True):
        return (3, layer_radii_combinations[3])

    elif (orfs == True and gc == True and custom_histogram == False):
        return (4, layer_radii_combinations[4])

    elif (orfs == True and gc == False and custom_histogram == True):
        return (5, layer_radii_combinations[5])

    elif (orfs == False and gc == True and custom_histogram == True):
        return (6, layer_radii_combinations[6])

    elif (orfs == True and gc == True and custom_histogram == True):
        return (7, layer_radii_combinations[7])



def create_circos_conf_file(output, karyotype_file, links_file, orfs, gc, custom_histogram, orf_forward_file=None, orf_reverse_file=None, gc_file=None, custom_histogram_file=None):
    """
    Creates the conf file for a circos diagram based on the type and amount of layers.

    Args:
        output (str): Output file path for the circos diagram (e.g. path/circos.conf)
        karyotype_file (str): File path to karyotype.txt file
        orfs (bool): An orf layer is inlcuded in the diagram (True or False)
        gc (bool): A GC layer is inlcuded in the diagram (True or False)
        custom_histogram (bool): A custom histogram layer is inlcuded in the diagram (True or False)
        orf_forward_file (str): File path to ORF.txt file
        orf_reverse_file (str): File path to ORF_reverse.txt file
        gc_file (str): File path to gc_content.txt if GC track supposed to be inlcuded\
        custom_histogram_file (str): File path to txt file data for the custom histogram layer
        links_file (str): File path to links.txt file

    Returns:
        None

    Notes:
        Layer config variable holds the layer combination, which are outlined below:
            0 = karyotype, links
            1 = karyotype, ORFs, links
            2 = karyotype, GC, links
            3 = karyotype, custom_histogram, links
            4 = karyotype, ORFs, GC, links
            5 = karyotype, ORFs, custom_histogram, links
            6 = karyotype, GC, custom_histogram, links
            7 = karyotype, ORFs, GC, custom_histogram, links
    """

    # Determine the radii for the given layer configuration
    layer_radii_and_type = find_layer_radii(orfs, gc, custom_histogram)
    layer_config = layer_radii_and_type[0]
    layer_radii_dict = layer_radii_and_type[1]

    # Erase everything in the conf file
    try:
        open(output, 'w').close()
    except Exception as e:
        print(e)
        print("Was not able to create the circos.conf configuration file")

    # Open the file in append mode
    try:
        f = open(output, "a")
    except Exception as e:
        print(e)
        print("Was not able to create the circos.conf configuration file.")


    f.write(
"""
################################################################
# Ideogram and Karyotype Setup
################################################################

karyotype = """ + karyotype_file + """

<ideogram>

<spacing>
default = 0.01r
</spacing>

# Ideogram position, fill and outline
radius    = 0.80r
thickness = 90p
fill      = yes
stroke_color = black
stroke_thickness = 8p

show_bands = yes
fill_bands = yes
band_transparency = 0
show_label = yes
# see etc/fonts.conf for list of font names
label_font = bold
label_radius = 1r + 130p
label_size = 55
label_parallel = yes
</ideogram>


################################################################
# Ticks and labels on fosmid karyotype bars
################################################################

show_ticks = yes
show_tick_labels = yes

<ticks>
skip_first_label = no
skip_last_label = no
radius = dims(ideogram,radius_outer)
# Setting multiplier to this means 10 000 will show as 10
multiplier = 1e-3
color = black
thickness = 4p
size = 20p

# 5kb bars with no label
<tick>
spacing = 5000u
size = 0.2r
show_label = no
label_size = 1r
thickness = 4p
color = black
</tick>

# 10kb bars with labels
<tick>
spacing = 10000u
size = 0.3r
show_label = yes
label_size = 1r
label_offset = 0.5r
thickness = 6p
color = black
</tick>
</ticks>

################################################################
# Beginning of Plots Section
################################################################
<plots>
"""
    )

    if layer_config in (1,4,5,7):
        f.write(
"""
################################################################
# ORF Forward and Reverse Layers
################################################################

# Forward Strand ORFs
<plot>
type            = tile
layers_overflow = hide
file        = """ + orf_forward_file + """
r1          = """ + layer_radii_dict['ORF'][1] + """
r0          = """ + layer_radii_dict['ORF'][0] + """
orientation = in
layers      = 3
margin      = 0.02u
thickness   = 15
padding     = 8
stroke_thickness = 1
stroke_color   = black
color = outer_orf_blue_color
<backgrounds>
<background>
color = outer_orf_background_color
</background>
</backgrounds>
</plot>

# Reverse Strand ORFs
<plot>
type            = tile
layers_overflow = hide
file        = """ + orf_reverse_file + """
r1          = """ + layer_radii_dict['ORF_reverse'][1] + """
r0          = """ + layer_radii_dict['ORF_reverse'][0] + """
orientation = in
layers      = 3
margin      = 0.02u
thickness   = 15
padding     = 8
stroke_thickness = 1
stroke_color     = black
color = inner_orf_red_color
<backgrounds>
<background>
color = inner_orf_background_color
</background>
</backgrounds>
</plot>
"""
        )
    if layer_config in (2,4,6,7):
        f.write(
"""
################################################################
# GC Content Layer
################################################################

<plot>
#Histogram Data Format: chr start end value [options]
type      = histogram
file      = """ + gc_file + """
r1        = """ + layer_radii_dict['GC'][1] + """
r0        = """ + layer_radii_dict['GC'][0] + """
max       = 100
min       = 0
stroke_type = outline
thickness   = 1
color       = vdgrey
extend_bin  = no
<backgrounds>
<background>
color = gc_histogram_background_color
</background>
</backgrounds>
<axes>
<axis>
spacing   = 0.1r
color     = vlgrey
thickness = 1
</axis>
</axes>
fill_color = gc_histogram_color
</plot>
"""
        )

    if layer_config in (3,5,6,7):
        f.write(
"""
################################################################
# Custom Histogram Layer
################################################################

<plot>
#Histogram Data Format: chr start end value [options]
type      = histogram
file      = """ + custom_histogram_file + """
r1        = """ + layer_radii_dict['custom_histogram'][1] + """
r0        = """ + layer_radii_dict['custom_histogram'][0] + """
max       = 100
min       = 0
stroke_type = outline
thickness   = 1
color       = vdgrey
extend_bin  = no
<backgrounds>
<background>
color = gc_histogram_background_color
</background>
</backgrounds>
<axes>
<axis>
spacing   = 0.1r
color     = vlgrey
thickness = 1
</axis>
</axes>
fill_color = gc_histogram_color
</plot>

"""
        )

    f.write(
"""
################################################################
# End of Plots Section
################################################################
</plots>

################################################################
# Links
################################################################
<links>
<link>
file = """ + links_file + """
radius = """ + layer_radii_dict['links'] + """
bezier_radius = 0.05r
thickness = 1
ribbon = yes
stroke_color = black_a4
stroke_thickness = 4
</link>
</links>

################################################################
# The remaining content is standard and required. It is imported
# from default files in the Circos distribution.
################################################################

<image>
<<include etc/image.conf>> # Included from Circos distribution.
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O and other system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>

################################################################
# Defining custom colors
################################################################
<colors>
outer_orf_background_color = 120,120,120,0.30
inner_orf_background_color = 120,120,120,0.2
gc_histogram_background_color = 120,120,120,0.1
outer_orf_blue_color = vdgrey
inner_orf_red_color = vdgrey
gc_histogram_color = vdgrey
</colors>
"""
    )
    f.close()


def prtn_domain_legend(karyotype, file_name, ncol):
    """
    Take a karyotype.txt file and produce a legend for all the unique bands

    Args:
        karyotype (str): File path to a karyotype.txt file
        file_name (str): File name to be saved (should be name.png)
        ncol (int): Number of coloumns of labels

    Returns:
        None
    """
    with open(karyotype, 'r') as f:

        domains = {}

        for line in f:
            if line.startswith("band"):
                domains[line.split()[2]] = line.split()[6].replace("rgb", "")

        domains_list = list(domains.keys())
        rgb = [tuple(map(int, color[1:-1].split(','))) for color in list(domains.values())]

        draw_prtn_domain_legend(domains_list, rgb, file_name, ncol)


def draw_prtn_domain_legend(labels, colors, file_name, ncol):
    """
    Takes in labels and corresponding rgb colors and produces a legend

    Args:
        labels (list<str>): A list of labels
        Colors (list<tuples>): A list of rbg tubles (e.g. [(255,10,11), (56,110,255)]
        file_name (str): File name to be saved (should be name.png)
        ncol (int): Number of coloumns of labels

    Returns:
        None
    """

    # Convert colors to 0-1 scale
    normalized_colors = []
    for color in colors:
        normalized_colors.append((color[0]/255, color[1]/255, color[2]/255))

    plt.figure(figsize=(12,6))
    f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
    plt.axis('off')
    handles = [f("s", normalized_colors[i]) for i in range(len(labels))]
    legend = plt.legend(handles, labels, loc=10, framealpha=1, frameon=False, ncol=ncol, mode="None", prop={'size': 6})

    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(file_name, dpi=400, bbox_inches=bbox)


def karyotype_legend(karyotype_file, output_file):
    """
    Takes in a list of fosmid names with their corresponding number (e.g. 1 : FosmidA7)
    and produces a legend.

    Args:
        karyotype_file (str): File path to karyotype.txt file containing fosmid information
        output_file (str): File path to file where the legend is to be saved to (e.g. <path>/name.png)

    Returns:
        None
    """

    fosmid_names = []

    with open(karyotype_file, 'r') as f:

        for line in f:
            if line.startswith("chr"):
                fosmid_names.append(line.split()[3] + ' : ' + line.split()[2])

    plt.figure(figsize=(12,6))
    f = lambda : plt.plot([],[],marker=None, color='white', ls="none")[0]
    plt.axis('off')
    handles = [f() for i in range(len(fosmid_names))]
    legend = plt.legend(handles, fosmid_names, loc=9, framealpha=1, frameon=False, mode="None", prop={'size': 6}, handletextpad=0, borderpad=1)

    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(output_file, dpi=400, bbox_inches=bbox)
