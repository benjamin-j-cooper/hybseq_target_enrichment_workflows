"""
This script take folder and find trees ending in users choice and trims non-overlapping taxa between those trees using Phix

Phix (pxrmt) must be on path

"""

import os, sys, glob, shutil
from Bio import Phylo
import pandas as pd


def remove_tips(DIR, tree_file_ending_1, tree_file_ending_2):
	""""Takes two trees and and trims non-overlapping taxa"""

	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"
	
	filecount = 0


	for i in os.listdir(DIR):
		if i.endswith(tree_file_ending_1):
			print i
			filecount += 1
			original_tre = Phylo.read(DIR+i, "newick")
			original_tre_tips = []
			for leaf in original_tre.get_terminals(): original_tre_tips.append(leaf.name)
			ori_tip_db = pd.DataFrame(original_tre_tips)
			
			
			path_tree, files_tree = os.path.split(i)
			tree_name = str(files_tree)
			base_name_tree = tree_name.split( "." )
			

			
			#if fail remove "." after base_name_tree[0]
			trim_tree = glob.glob(DIR + base_name_tree[0]+ "." + '*.'+ tree_file_ending_2)
			
			print trim_tree
			
			trimmed_tre = Phylo.read(trim_tree[0], "newick")
			trimmed_tre_tips = []
			for leaf in trimmed_tre.get_terminals(): trimmed_tre_tips.append(leaf.name)
			trim_tip_db = pd.DataFrame(trimmed_tre_tips)
			
			
			tips_concat = pd.concat([ori_tip_db, trim_tip_db])
			tips_concat.columns = ['tip_labels']
			uniques = tips_concat.drop_duplicates(keep=False)
			list = uniques['tip_labels'].values.tolist()
			
			if len(list) == 0:
				shutil.copy2(DIR+i, (DIR+base_name_tree[0]+"."+base_name_tree[1]+".pxrmt.tm"))
			
			else:
			
				cmd = ["pxrmt","-t", DIR+i, "-o", DIR+base_name_tree[0]+"."+base_name_tree[1]+".pxrmt.tm", "-n"]
				print ((" ".join(cmd))+(" ")+(','.join(list)))
				os.system((" ".join(cmd))+(" ")+(','.join(list)))
            
            
	assert filecount > 0, \
		"No file end with "+tree_file_ending+" found in "+DIR


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python diff_two_trees.py DIR tree_file_ending_1 tree_file_ending_2 "
		sys.exit(0)

	DIR,tree_file_ending_1,tree_file_ending_2 = sys.argv[1:]
	remove_tips(DIR,tree_file_ending_1, tree_file_ending_2)
