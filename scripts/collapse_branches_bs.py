###This script takes a folder of trees (newick) and collapse branches given a specific bootstrap threshold


import os,sys
from Bio import Phylo


def collapse_trees(inDIR,file_ending,bs_min_value,outDIR):

	if os.path.isabs(inDIR) == False: inDIR = os.path.abspath(inDIR)
	if inDIR[-1] != "/": inDIR += "/"
	if os.path.isabs(outDIR) == False: outDIR = os.path.abspath(outDIR)
	if outDIR[-1] != "/": outDIR += "/"
	if outDIR == ".": outDIR = os.getcwd()
	
	
	filecount = 0


	for i in os.listdir(inDIR):
		if i.endswith(file_ending):
			print i
			filecount += 1

			tree = Phylo.read(inDIR+i, 'newick')
			tree.collapse_all(lambda c: c.confidence is not None and c.confidence < int(bs_min_value))
			basename = os.path.splitext(i)[0]
			new_tree= outDIR+basename+".col_"+bs_min_value+".tre"
			Phylo.write(tree, new_tree, 'newick')
			
	assert filecount > 0, \
		"No file end with "+file_ending+" found in "+inDIR


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "Usage:"
		print "python collapse_branches_bs.py inDIR, file_ending ,bs_min_value, outDIR"
	else:	
		collapse_trees(inDIR=sys.argv[1],file_ending=sys.argv[2],bs_min_value=sys.argv[3],outDIR=sys.argv[4])
	sys.exit(0)		