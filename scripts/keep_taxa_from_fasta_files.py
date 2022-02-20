"""
This script take a folder of fasta files and 
write fasta files with only the taxa on a given list

"""

import os, sys, re
from Bio import SeqIO

def write_fasta(fasta_DIR, fasta_file_ending, taxa_to_keep_file, outDIR):
	""""Takes a fasta file and tree and write fasta file only the taxa from the tree"""

	if os.path.isabs(fasta_DIR) == False: fasta_DIR = os.path.abspath(fasta_DIR)
	if fasta_DIR[-1] != "/": fasta_DIR += "/"
	
	if os.path.isabs(outDIR) == False: outDIR = os.path.abspath(outDIR)
	if outDIR[-1] != "/": outDIR += "/"
	if outDIR == ".": outDIR = os.getcwd()
	
	filecount = 0

	with open(taxa_to_keep_file) as f:
		taxa_to_keep = f.read().splitlines()
		#print taxa_to_keep
	

	    
	for i in os.listdir(fasta_DIR):
		if i.endswith(fasta_file_ending):
			print i
			filecount += 1
			original_fasta = SeqIO.parse(fasta_DIR+i, "fasta")
			original_fasta_names = []
			for record in original_fasta: original_fasta_names.append(record.id)
			#print original_fasta_names
			
			path_fasta, files_fasta = os.path.split(i)
			fasta_name = str(files_fasta)
			base_name_fasta = fasta_name.split( "." )
			#print base_name_fasta
			
			filtered_fasta_names = []
	        for taxa in taxa_to_keep:
	            r = re.compile(str(taxa)+"@.*")
	            filtered_fasta_names.append(filter(r.match, original_fasta_names))
	            
	        #print filtered_fasta_names
	        
	        flatlist = []
	        for sublist in filtered_fasta_names:
	        	for item in sublist:
	        		flatlist.append(item)
	        	#print flatlist

	
			#if len(flatlist) == len(taxa_to_keep):
			if len(flatlist) >=3:

				new_fasta_file = outDIR+base_name_fasta[0]+".ortho.rd.fa"
				#print new_fasta_file
				file=open(new_fasta_file,'w')
	
				for record in SeqIO.parse(fasta_DIR+i,"fasta"):
				    if record.id in flatlist:
				        id =record.id
				        seq=str(record.seq)
				        file.write(">" + id + "\n" + seq + "\n")
				file.close()
				

	assert filecount > 0, \
		"No file end with "+fasta_file_ending+" found in "+fasta_DIR



if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python keep_taxa_from_fasta_files_complete_taxa.py inDIR fasta_file_ending taxa_to_keep_txt_file outDIR"
		sys.exit(0)

	fasta_DIR, fasta_file_ending, taxa_to_keep_file, outDIR  = sys.argv[1:]
	write_fasta(fasta_DIR, fasta_file_ending, taxa_to_keep_file, outDIR)
