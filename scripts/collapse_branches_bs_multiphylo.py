###This script takes a files of trees and collapse branches given a specific bootstrap threshold


import os,sys,io
from Bio import Phylo


def collapse_trees(inMultiTree,bs_min_value):

	
	basename = os.path.splitext(inMultiTree)[0]
	cwd = os.getcwd()
	new_mulititree= cwd+"/"+basename+".col_"+bs_min_value+".tre"		
	f= open(new_mulititree,"w")
	
	trees = Phylo.parse(inMultiTree, 'newick')
	
	for tree in trees:
		#print("Total branch length %0.2f" % tree.total_branch_length())
		tree.collapse_all(lambda c: c.confidence is not None and c.confidence < int(bs_min_value))
		#print("Total branch length %0.2f" % tree.total_branch_length())
		Phylo.write(tree, f, 'newick')
	
	f.close()	
			
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "Usage:"
		print "python collapse_branches_bs.py inMultiTree, bs_min_value"
	else:	
		collapse_trees(inMultiTree=sys.argv[1],bs_min_value=sys.argv[2])
	sys.exit(0)		