import fosvis

contigs = '/mnt/nfs/sharknado/Sandbox/Tylo/circos_outputs_and_sample_data/sample_data/CD1_contigs.fasta'
output_dir = '/mnt/nfs/sharknado/Sandbox/Tylo/circos_outputs_and_sample_data'
hmm_db = '/mnt/nfs/sharknado/Sandbox/Tylo/databases/hmmer/hmmscan/Pfam-A.hmm'
proj_name = 'CD1_contigs_3'

print("Running...")
fosvis.create_circos_data(contigs, output_dir, proj_name, hmm_db, e_value=0.0001, min_contig_length=10000, min_blast_similarity_length=400, min_blast_similarity_percentage=90, include_domains=False, gc=True, gc_interval_len=200, blast_type='blastn', keep_blast_data=False)
fosvis.make_diagram(output_dir + "/" + proj_name + "/circos_diagram_input_data", 5)
print("Finished Running...")
