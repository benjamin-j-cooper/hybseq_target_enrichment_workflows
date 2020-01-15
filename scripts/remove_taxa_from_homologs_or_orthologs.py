import os, sys, re, shutil
from Bio import Phylo


def remove_taxa(inDIR, tree_file_ending, outDIR, taxa_to_remove_file):
	""""Takes a list of taxa and remove it from trees"""

	if os.path.isabs(inDIR) == False: inDIR = os.path.abspath(inDIR)
	if inDIR[-1] != "/": inDIR += "/"
	
	if os.path.isabs(outDIR) == False: outDIR = os.path.abspath(outDIR)
	if outDIR[-1] != "/": outDIR += "/"
	
	if outDIR == ".": outDIR = os.getcwd()

	filecount = 0
	
	with open(taxa_to_remove_file) as f:
	    taxa_to_remove = f.read().splitlines()
	    print taxa_to_remove

	for i in os.listdir(inDIR):
	    if i.endswith(tree_file_ending):
	    	print i
	       
	    	filecount += 1
	    	original_tre = Phylo.read(inDIR+i, "newick")
	    	original_tre_tips = []
	    	for leaf in original_tre.get_terminals(): original_tre_tips.append(leaf.name)
	        
	    	path_tree, files_tree = os.path.split(i)
	    	tree_name = str(files_tree)
	    	base_name_tree = tree_name.split( "." )
	    	
	    	newlist = []
	    	for taxa in taxa_to_remove:
	    		r = re.compile(str(taxa)+".*")
	    		newlist.append(filter(r.match, original_tre_tips))
	    		
	    	flatlist = []
	    	for sublist in newlist:
	    		for item in sublist:
	    			flatlist.append(item)

            #print flatlist
            
            if len(flatlist) == 0:
            	shutil.copy2(inDIR+i, (outDIR+base_name_tree[0]+"."+base_name_tree[1]+".txrm.tm"))
                
            else:
                cmd = ["pxrmt","-t", inDIR+i, "-o", outDIR+base_name_tree[0]+"."+base_name_tree[1]+".txrm.tm", "-n"]
                print ((" ".join(cmd))+(" ")+(','.join(flatlist)))
                os.system((" ".join(cmd))+(" ")+(','.join(flatlist)))
            
            
	assert filecount > 0, \
		"No file end with "+tree_file_ending+" found in "+inDIR

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python remove_taxa_from_homologs_or_orthologs.py inDIR tree_file_ending outDIR taxa_to_remove_file"
		sys.exit(0)

	inDIR, tree_file_ending, outDIR, taxa_to_remove_file = sys.argv[1:]
	remove_taxa(inDIR, tree_file_ending, outDIR, taxa_to_remove_file)




