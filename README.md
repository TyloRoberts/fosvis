Author: Tylo Roberts  
Contact: tylojroberts@gmail.com

# FosVis Overview  
A pip-installable python package designed for visualization of fosmids from functional screens.  Inputs are contigs from an assembly (representing fosmids) and the output is a circos diagram.  The diagram shows sequence or protein homology, open reading frames, GC content and protein domains located on the fosmids.  The package is designed for the following workflow:

![alt text](https://github.com/TyloRoberts/fosvis/blob/master/images/fosvis_context.png?raw=true)

# Requirements
* Python3 (must be run with python3)
* Pip-installable packages: BioPython, Pandas, Numpy, Matplotlib, seaborn, tqdm
* External Software (all need to be appended to PATH)
	* circos (v0.69-8)
	* prodigal (v2.6.3)
	* blastn (v2.9.0)
	* bl2seq (includes tblastx)
	* hmmscan (HMMER 3.3)

* To install circos using conda: ```conda install circos```
	* You will need to have already set up the bioconda channel: ```conda config --add channels bioconda```
	* See more at <https://bioconda.github.io/recipes/circos/README.html>


# Getting Started
The package has two main methods that control the functionality.  To create the data for the image run ```create_circos_data()```.  Then run ```make_diagram()``` using the created data to output the circos diagram.  The following provides the most basic use of the package:

1. Ensure all requirements are met
2. ```pip install fosvis```
3. Run the following example script with ```python3```

```
from fosvis import fosvis

# fasta contigs file (all contigs in one file)
contigs = '<path_to_contigs_fasta_file>'

# File path to a directory where the project folder will be created
output_dir = '<path_to_an_output_directory>'

# A Project Name (output directory will be created with this name)
proj_name = 'test_project'

print("Running...")
fosvis.create_circos_data(contigs, output_dir, proj_name)
fosvis.make_diagram(output_dir + "/" + proj_name + "/circos_diagram_input_data")
print("Finished Running...")
```

You can now find your output image in the following directory: ```output_dir/proj_name/circos_diagram_input_data/circos_diagram.png```

# Additional Parameters

### ```create_circos_data()```

```
def create_circos_data(contigs, output_dir, project_title, hmm_db='', e_value=0.01, min_contig_length=10000,
 min_blast_similarity_percentage=90, min_blast_similarity_length=300, link_transperacny='0.60',
 percent_link_overlap_tolerance=50, include_domains=False, gc=True, gc_interval_len=100, blast_type='blastn', keep_blast_data=False)
```

**Required Parameters**

* ```contigs``` (str): File path to a FASTA file containing contigs
* ```output_dir``` (str): File path to a directory where the project folder will be created
* ```project_title``` (str): Name of the project directory that is created in the ouput_dir (and used for some file prefixes)


**Parameters with Default Values**

* ```hmm_db``` (str): File path to a pressed hmm database (default='')
* ```e_value``` (int): The e-value used for the hmmscan search (default=0.01)
* ```min_contig_length``` (int): The minimum size of a contig to be used in the diagram (idea is that each contig should represent one fosmid) (default=10000)
* ```min_blast_similarity_percentage``` (int): Min nucleotide percent similarity (blast percent identity) for link data (default=90)
* ```min_blast_similarity_length``` (int): Min nucleotide similarity length (blast align length) for link data (default=300)
* ```link_transperacny``` (**str**): The transparency of the links in the diagram in range 0 - 1 (0 is transparent) (default='0.60')
* ```percent_link_overlap_tolerance``` (int): The amount of overlap percentage required to classify two links as similar enough to color them with the same color (default=50)
* ```inlcude_domains``` (boolean): If True will create protein domain band data, if False will not (default=False)
	* If set to True, you must also set the ```hmmdb``` parameter to a valid hmm database
* ```gc``` (boolean): If True, the diagram will include a track for GC content, if False it will not (default=True)
* ```gc_interval_len``` (int): The interval over which gc content is calculated for the histogram (default=100)
* ```blast_type``` (str): The type of blast to use to create the links (can be 'blastn' OR 'tblastx') (default='blastn')
* ```keep_blast_data``` (bool): If True will keep the raw blast data, won't if False (default=False)


### ```make_diagram()```
```
make_diagram(data_dir, ncol=2)
```
**Required Parameters**

* ```data_dir``` (str): File path to a directory containing the following files created by ```create_circos_data()```:
	* ORF.txt  
	* ORF_reverse.txt  
	* Links.txt  
	* Karyotype.txt  
	* circos.conf  

* i.e. it is the path to the ```circos_diagram_input_data``` directory created within the directory created by ```create_circos_data()```


**Parameters with Default Values**

* ```ncol``` (int): Number of columns in the protein domain legend if applicable (default=2)

# Output Visualization
![alt text](https://github.com/TyloRoberts/fosvis/blob/master/images/fosvis_sample_image.png?raw=true)

The circos diagram output image has the following main features (depending on the parameters used):

* Outer layer: Represents length/position of fosmid sequences (scale bars show fosmid size in kb)
* Coloured Bands on Outer Layer: Represent protein domains with each unique domain having an associated color (shown in a separate legend) (will only be shown if ```inlcude_domains=True```)
* 2 Tile Layers (layers 2 & 3): Represent open reading frames (2nd layer shows forward strand, 3rd layer shows reverse strand)
* Histogram Layer: GC content (only included if ```gc=True```)
* Inner Ribbons: The ribbons (links) represent nucleotide/protein homology using blastn or tblastx between fosmids
	* Links representing the same (or part of the same) sequence are grouped by color


# Output Files
The script outputs the following directories and files after running ```create_circos_data()``` and ```make_diagram()```:

**```intermediate_ouputs/prodigal``` Directory**

* ```prodigal_prtn_seq_output.faa``` - Protein sequences of prodigal open reading frames identified in the fosmids
* ```prodigal_orf_output```- Prodigal open reading frame position data

**```intermediate_ouputs/hmmscan``` Directory**

* Empty unless ```include_domains=True```
* ```<project_title>_tblout.txt```  - hmmscan tbl format output file
* ```<project_title>_out.txt``` - hmmscan full format output file
* ```<project_title>_domtblout.txt``` - hmmscan domtbl format output file

**```intermediate_ouputs/correct_length_contigs``` Directory**

* ```<contig_name>.fasta``` - An individual fasta file for each contig in the input contigs file that is >= to the minimum contig size provided
* ```all_correct_size_contigs.fasta``` - All correct sized contigs in one fasta file

**```intermediate_ouputs/blastn``` Directory**

* Contains blastn outputs if ```keep_blast_data=False```

**```circos_diagram_input_data``` Directory**

* ```ORF.txt``` - Circos open reading frame (forward strand) input data
* ```ORF_reverse.txt``` - Circos open reading frame (reverse strand) input data
* ```links.txt``` - Circos sequence similarity link data
* ```karyotype.txt``` - Circos outer layer input data
* ```protein_domain_legend.png``` - Legend showing protein domain colors shown on the circos diagram (Empty unless ```include_domains=True```)
* ```karyotype_legend.png``` - Legend for fosmid names represented as 1,2,3 etc. on the diagram
* ```circos.conf``` - Configuration file for circos software
* ```circos_diagram.svg``` - Circos diagram output (svg format)
* ```circos_diagram.png``` - Circos diagram output (png format)


# Manual Image Modifications
Manual modifications of the data are sometimes necessary to clean up the image.  In order to manually edit the image, simply manually edit the circos input files and re-run ```make_diagram()```. A general outline of how to make all modifications is described first and specifics on common modifications are described after.

**How to Make Modifications**  
1. Open one of the 4 circos input txt files with a text editor  
2. Make the modification (details of common modifications described below)  
3. Save the changes to the same file  
4. Run the ```make_diagram()``` function with the ```circos_diagram_input_data``` directory (containing the modified file)  
5. The ```circos_diagram.png``` file and legends will be overwritten with the changes made to the input files  

**Details on Common Modifications (What to do in step 2 above)**  

1. Removing irrelevant protein domain annotations
	* Initially, the protein domain legend is likely filled with irrelevant protein domains that need to be removed
	* In the karyotype.txt file the protein domain bands are represented by lines that start with ‘band’
	* In karyotype.txt the second column is the fosmid the band is on, the third is a brief name of the domain and the fourth column is a more verbose name of the protein domain (Note: spaces from database name are replaced with underscores)
	* For any protein domain that is irrelevant to your diagram, delete the line and when make_diagram() is run again the circos diagram and legend will reflect the changes

2. Removal of links and ORFs
	* Open the links.txt or ORF.txt file and simply delete the row containing the link or ORF that you want to remove

3. Changing colors of links (same for protein domains)
	* In the links.txt file the color of each link is denoted in the final column after the ‘color=’ tag and is surrounded in brackets
		* Color is denoted in the format (r,g,b,luminosity)
	* Modify the rgb/luminosity values for each link that you would like to change the color of
		* Keep in mind that some links are already grouped so be sure to change all the colors in the group to maintain those groups
		* If you search the file for the color that is currently entered it will show you all the occurrences of that color and hence all the links that are grouped with that link
	* Note about protein domains:
		* Same but in format ‘rgb(r,g,b)’ and have no luminosity
		* Recurrent domains already have same color so if you change one change all to that color otherwise the legend won’t be accurate

4. Changing Fosmid Labels from 1,2,3 etc.
	* Open the karyotype.txt file and the chromosomes are represented by lines starting with ‘chr’
	* The label for the fosmid shown on the diagram is in the fourth column (including the ‘-’ as a column)
	* Change that value to the desired label and re-run make_diagram()

5. Changing Diagram Layout/Features (not recommended)
	* To change the fundamental layout of the image (layer size/positions, adding additional layers) you need to edit the circos.conf file
	* See the circos documentation for more information: <http://circos.ca/documentation/>

# Implementation Details
The script uses a variety of tools to create the various data inputs that the circos software can use to create an image.  The main tools are shown in the diagram below and some additional details are provided about notable configurations of those tools.
![alt text](https://github.com/TyloRoberts/fosvis/blob/master/images/implementation.png?raw=true)



**Links**

* Every combination of pairs of fosmids are blasted against each other (using blastn or tblastx)
* The outputs are parsed and only blast alignments are only kept if the sequences have a percent similarity >= the parameter min_blast_similarity_percentage and a length >= the parameter ```min_blast_similarity_length``` given in ```create_circos_data()```
* A function is then run that looks for sets of links that are potentially representing the same (or some of the same) sequence and groups them by assigning them all with the same color

**prodigal**

* Some prodigal settings worth noting are:
	* The -c flag was not used so genes were allowed to run off the edge of the fosmid
	* Prodigal also uses GTG as a start codon


# Small Details Worth Noting

* For the ORF layer, if more than 3 overlap, then any more are not shown in the image
* Whether the link pinches in the middle or not is as a result of link position data having the orign\_start > or < the origin\_end and the terminus\_start > or < the terminus\_end
* The chromosome name used will be the sequence id from the fasta file up to (not including) the first space if there is a space in the fasta sequence id.  Thus, make sure for every sequence every header is unique up to the first space.
* Circos indexing:
	* The circos software takes each bp position to be a range
	* If your karyotype starts at 1 (as in fosvis) the first position would be the range from 1 - 2, the 2nd would be the range from 2 - 3


## References

**Software**

*Hyatt, D., Chen, G., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119*

*Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990;215(3):403-410. doi:10.1016/S0022-2836(05)80360-2*

*HMMER. (2020). http://hmmer.org/.*

*Krzywinski, M., Schein, J., Birol, I., Connors, J., Gascoyne, R., & Horsman, D. et al. (2009). Circos: An information aesthetic for comparative genomics. Genome Research, 19(9), 1639-1645. doi: 10.1101/gr.092759.109*



**Sample Image Data**

*Mewis K, Armstrong Z, Song YC, Baldwin SA, Withers SG, Hallam SJ. Biomining active cellulases from a mining bioremediation system. J Biotechnol. 2013;167(4):462-471. doi:10.1016/j.jbiotec.2013.07.015*

## Acknowledgements

UBC Hallam Lab  
Avery Noonan  
Connor Morgan-Lang
