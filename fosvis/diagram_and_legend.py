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



def create_circos_conf_file(output, karyotype, links, ORF, ORF_reverse, gc_data=None):
    """
    Creates the conf file for a circos diagram.

    Args:
        output (str): Output file path for the circos diagram (e.g. path/circos.conf)
        karyotype (str): File path to karyotype.txt file
        links (str): File path to links.txt file
        ORF (str): File path to ORF.txt file
        ORF_reverse (str): File path to ORF_reverse.txt file
        gc_data (str): File path to gc_content.txt if GC track supposed to be inlcuded

    Returns:
        None
    """

    try:
        f = open(output, "w")
    except Exception as e:
        print(e)
        print("Was not able to create the circos.conf configuration file.")

    if gc_data == None:
        f.write(
"""# circos.conf configuration file

karyotype = """ + karyotype + """

<ideogram>

<spacing>
default = 0.01r
</spacing>

radius    = 0.85r
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

<tick>
spacing = 5000u
size = 0.2r
show_label = no
label_size = 1r
thickness = 4p
color = black
</tick>

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

<plots>
type            = tile
layers_overflow = hide

<plot>
file        = """ + ORF + """
r1          = 0.96r
r0          = 0.91r
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

<plot>
file        = """ + ORF_reverse + """
r1          = 0.88r
r0          = 0.83r
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

</plots>

<links>
<link>
file = """ + links + """
radius = 0.79r
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

<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Custom defined colors
<colors>

outer_orf_background_color = 120,120,120,0.30
inner_orf_background_color = 120,120,120,0.2
outer_orf_blue_color = vdgrey
inner_orf_red_color = vdgrey

</colors>

# Debugging, I/O and other system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>"""
            )
        f.close()
    # Inlcuding GC layer
    else:
        f.write(
"""# circos.conf configuration file - gc inlcuded

karyotype = """ + karyotype + """

<ideogram>

<spacing>
default = 0.01r
</spacing>

radius    = 0.85r
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

<tick>
spacing = 5000u
size = 0.2r
show_label = no
label_size = 1r
thickness = 4p
color = black
</tick>

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

<plots>
type            = tile
layers_overflow = hide

<plot>
file        = """ + ORF + """
r1          = 0.96r
r0          = 0.91r
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

<plot>
file        = """ + ORF_reverse + """
r1          = 0.88r
r0          = 0.83r
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

<plot>

#Histogram Data Format: chr start end value [options]

type      = histogram
file      = gc_content.txt

r1        = 0.8r
r0        = 0.7r
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

</plots>

<links>
<link>
file = """ + links + """
radius = 0.66r
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

<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Custom defined colors
<colors>

outer_orf_background_color = 120,120,120,0.30
inner_orf_background_color = 120,120,120,0.2
gc_histogram_background_color = 120,120,120,0.1
outer_orf_blue_color = vdgrey
inner_orf_red_color = vdgrey
gc_histogram_color = vdgrey

</colors>

# Debugging, I/O and other system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>"""
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
