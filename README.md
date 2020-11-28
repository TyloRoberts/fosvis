Author: Tylo Roberts  
Contact: tylojroberts@gmail.com

##Overview  
A python module designed for visualization of fosmids from functional screens.  Inputs are contigs from an assembly (representing fosmids) and the output is a circos diagram.  The diagram shows sequence or protein homology, open reading frames, GC content and protein domains located on the fosmids.

##Circos Output Image
The circos diagram output image has the following main features:

* Outer layer: Represents the fosmid sequences (scale bars show fosmid size in kb)
* Coloured Bands on Outer Layer: Represent protein domains with each unique domain having an associated color (shown in a separate legend)
* Layer 2 and 3: Represent open reading frames (2nd layer shows forward strand, 3rd layer shows reverse strand)
* Layer 4: GC content
* Inner Ribbons: The ribbons (links) represent nucleotide/protein 	homology using blastn or tblastx between fosmids
	* Links representing the same (or part of the same) sequence are grouped by color

##Requirements
* Python3
* Pip packages: BioPython, Pandas, Numpy, Matplotlib, seaborn, tqdm
* External Software (all need to be appended to PATH)
	* circos (v0.69-8)
	* prodigal (v2.6.3)
	* blastn (v2.9.0)
	* hmmscan (HMMER 3.3)

* To install circos using conda: ```conda install circos```
	* You will need to have already set up the bioconda channel: ```conda config --add channels bioconda```

##Use
The module has two main methods that control the functionality.  To create the data for the image run ```create_circos_data()```.  Then run ```make_diagram()``` using the created data to output the circos diagram.  The two functions are described below:  


```
def create_circos_data(contigs, output_dir, project_title, hmm_db, e_value=0.01, min_contig_length=30000,
 min_blast_similarity_percentage=90, min_blast_similarity_length=300, link_transperacny='0.60',
 percent_link_overlap_tolerance=50, include_domains=True, gc=True, gc_interval_len=100, blast_type='blastn', keep_blast_data=False):
```

 <span style="color:blue"> split these up into required and default params</span>.

* contigs (str): File path to a FASTA file containing contigs
* output_dir (str): File path to a directory where the project folder will be created
* project_title (str): Name of the project directory that is created in the ouput_dir (and used for some file prefixes)
* hmm_db (str): File path to a pressed hmm database
* e_value (int): The e-value used for the hmmscan search
* min_contig_length (int): The minimum size of a contig to be used in the diagram (idea is that each contig should represent one fosmid)
* min_blast_similarity_percentage (int): Min nucleotide percent similarity (blast percent identity) for link data
* min_blast_similarity_length (int): Min nucleotide similarity length (blast align length) for link data
* link_transperacny (str): The transparency of the links in the diagram in range 0 - 1 (0 is transparent)
* percent_link_overlap_tolerance (int): The amount of overlap percentage required to classify two links as similar enough to color them with the same color
* inlcude_domains (boolean): If True will create protein domain band data, if false will not
* gc (boolean): If True the diagram will include a track for GC content, if False it will not
* gc_interval_len (int): The interval over which gc content is calculated for the histogram
* blast_type (str): The type of blast to use to create the links (can be 'blastn' OR 'tblastx')
* keep_blast_data (bool): If True will keep the raw blast data, won't if false

Note: If using inlcude_domains=False, can just set the hmmdb and e_value parameters to None - this doesnt work actually bc you concatenate it so just use ""

```
make_diagram(data_dir, ncol)
```
Description: Creates the circos diagram using the files in the data_dir

* data_dir (str): File path to a directory containing the following files:
	* ORF.txt  
	* ORF_reverse.txt  
	* Links.txt  
	* Karyotype.txt  
	* circos.conf  
* ncol (int): Number of columns in the protein domain legend  

##Example Use in a Python Script
```
from fosvis import fosvis

contigs = '<path>/<fosmids>.fasta'
output_dir = '<path>/<a_directory_to make_project_in>'
hmm_db = '<path>/Pfam-A.hmm'
proj_name = ‘test_project_name’

print(“Running…”)
fosvis.create_circos_data(contigs, output_dir, proj_name, hmm_db, e_value=0.0001, min_contig_length=1, min_blast_similarity_length=1, min_blast_similarity_percentage=1, include_domains=False, gc=True, gc_interval_len=500, blast_type='blastn', keep_blast_data=True)
fosvis.make_diagram(output_dir + "/" + proj_name + "/circos_diagram_input_data", 5)
print("Finished Running...")
```

##Outputs
The script outputs the following file structure for every project after running create_circos_data() and make_diagram():

 <span style="color:red"> image</span>

**prodigal Directory**

* prodigal\_prtn_seq_output.faa -> Protein sequences of prodigal open reading frames identified in the fosmids
* prodigal\_orf\_output -> Prodigal open reading frame position data

**hmmscan Directory**

* \<project\_title>_tblout.txt  -> hmmscan tbl format output file
* \<project\_title>_out.txt -> hmmscan full format output file
* \<project\_title>_domtblout.txt -> hmmscan domtbl format output file

**correct\_length\_contigs Directory**

* \<contig_name>.fasta -> An individual fasta file for each contig in the input contigs file that is >= to the minimum contig size provided
* all\_correct\_size\_contigs.fasta -> All correct sized contigs in one fasta file

**circos_diagram_input_data Directory**

* ORF.txt -> Circos open reading frame (forward strand) input data
* ORF_reverse.txt -> Circos open reading frame (reverse strand) input data
* links.txt -> Circos sequence similarity link data
* karyotype.txt -> Circos outer layer input data
* protein_domain_legend.png -> Legend showing protein domain colors shown on the circos diagram
* karyotype_legend.png -> Legend for fosmid names represented as 1,2,3 etc. on the diagram
* circos.conf -> Configuration file for circos software
* circos_diagram.svg -> Circos diagram output (svg format)
* circos_diagram.png -> Circos diagram output (png format)

Note: \<include note on blast file directory>

##Manual Image Modifications
Manual modifications of the data are sometimes necessary to clean up the image.  In order to manually edit the image, simply manually edit the circos input files and re-run make_diagram(). A general outline of how to make all modifications is described first and specifics on common modifications are described after.

**How to Make Modifications**  
1. Open the one of the 4 circos input txt files with a text editor  
2. Make the modification (details of common modifications described below)  
3. Save the changes to the same file  
4. Run the make_diagram() function with the circos_diagram_input_data directory (containing the modified file)  
5. The circos_diagram.png file and legends will be overwritten with the changes made to the input files  

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
	* See the circos website for more information

##Implementation Details
The script uses a variety of tools to create the various data inputs that the circos software can use to create an image.  The main tools are shown in the diagram below and some additional details are provided about notable configurations of those tools.

 <span style="color:red"> image</span>

**Links**

* Every combination of pairs of fosmids are blasted against each other (using blastn)
* The outputs are parsed and only blast alignments are only kept if the sequences have a percent similarity >= the parameter min_blast_similarity_percentage and a length >= the parameter min_blast_similarity_length given in create_circos_data()
* A function is then run that looks for sets of links that are potentially representing the same (or some of the same) sequence and groups them by assigning them all with the same color
* <add in tblastx info>

**prodigal**

* Some prodigal settings worth noting are:
	* The -c flag was not used so genes were allowed to run off the edge of the fosmid
	* Prodigal also uses GTG as a start codon


##Small Details Worth Noting

* For the ORF layer, if more than 3 overlap then any more are not shown in the image
* Whether the link pinches in the middle or not is as a result of link position data having the orign_start > or < the origin_end and the terminus_start > or < the terminus_end
* The chromosome name used will be the sequence id from the fasta file up to (not including) the first space if there is a space in the fasta sequence id.  Thus, make sure for every sequence every header is unique up to the first space.


 <span style="color:red"> refrences, circos indexing etx</span>
