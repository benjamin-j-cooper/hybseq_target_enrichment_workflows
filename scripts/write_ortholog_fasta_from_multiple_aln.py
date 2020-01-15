"""
This script take a folder of fasta files and 
write fasta files given a folder of trees 
and a minimum number of species
"""

import os, sys, glob
from Bio import Phylo, SeqIO

def write_fasta(fasta_DIR, tree_DIR, fasta_file_ending, tree_file_ending, outDIR):
	""""Takes a fasta file and tree and write fasta file only the taxa from the tree"""

	if os.path.isabs(tree_DIR) == False: tree_DIR = os.path.abspath(tree_DIR)
	if tree_DIR[-1] != "/": tree_DIR += "/"
	
	if os.path.isabs(fasta_DIR) == False: fasta_DIR = os.path.abspath(fasta_DIR)
	if fasta_DIR[-1] != "/": fasta_DIR += "/"
	
	if os.path.isabs(outDIR) == False: outDIR = os.path.abspath(outDIR)
	if outDIR[-1] != "/": outDIR += "/"
	if outDIR == ".": outDIR = os.getcwd()
	
	filecount = 0


	for i in os.listdir(tree_DIR):
		if i.endswith(tree_file_ending):
			print i
			filecount += 1
			original_tre = Phylo.read(tree_DIR+i, "newick")
			original_tre_tips = []
			for leaf in original_tre.get_terminals(): original_tre_tips.append(leaf.name)
			#print original_tre_tips
			
			path_tree, files_tree = os.path.split(i)
			tree_name = str(files_tree)
			base_name_tree = tree_name.split( "." )
			#print base_name_tree

			fasta_file = glob.glob(fasta_DIR + base_name_tree[0] + '*.'+ fasta_file_ending)
			#print fasta_file

			new_fasta_file = outDIR+base_name_tree[0]+".ortho.fa"
			#print new_fasta_file
			file=open(new_fasta_file,'w')

			for record in SeqIO.parse(fasta_file[0],"fasta"):
			    if record.id in original_tre_tips:
			        id =record.id
			        species_id,copy_id = id.split("@")
			        seq=str(record.seq)
			        file.write(">" + species_id + "\n" + seq + "\n")
			file.close()

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "python write_ortholog_fasta_from_multiple_aln.py fasta_DIR tree_DIR fasta_file_ending(no_dot) tree_file_ending(not_dot) outDIR"
		sys.exit(0)

	fasta_DIR, tree_DIR, fasta_file_ending, tree_file_ending, outDIR = sys.argv[1:]
	write_fasta(fasta_DIR, tree_DIR, fasta_file_ending, tree_file_ending, outDIR)
