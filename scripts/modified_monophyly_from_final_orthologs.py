from Bio import Phylo
import pandas as pd 
import mask_tips_by_taxonID_genomes
import newick3,phylo3,os,sys
from cStringIO import StringIO



def rename_tips(tree_DIR, tree_file_ending, clade_assignments, outDIR):
	"Replaces tip names in tree given a list of species with a clade assignment"

	if os.path.isabs(tree_DIR) == False: tree_DIR = os.path.abspath(tree_DIR)
	if tree_DIR[-1] != "/": tree_DIR += "/"
		
	if os.path.isabs(outDIR) == False: outDIR = os.path.abspath(outDIR)
	if outDIR[-1] != "/": outDIR += "/"
	if outDIR == ".": outDIR = os.getcwd()
	
	filecount = 0

	for x in os.listdir(tree_DIR):
		if x.endswith(tree_file_ending):
			print x
			filecount += 1

	
			df = pd.read_table(clade_assignments, delim_whitespace=True, header=None, names=['Clade','Species'])
		
			dict = {}
			for i in df['Species']:
				dict[i] = [df['Clade'][j] for j in df[df['Species']==i].index]

			path_tree, files_tree = os.path.split(x)
			tree_name = str(files_tree)
			base_name_tree = tree_name.split( "." )
		
			tree_DIR+x
			file = open(tree_DIR+x, "r")
			for line in file:
				tree = str(line)
				
			for key in dict.iterkeys():
				tree = (tree.replace(key, dict[key][0]))
			
			intree = mask_tips_by_taxonID_genomes.newick3.parse(tree)
			
				
			parsed = newick3.tostring(mask_tips_by_taxonID_genomes.monophyly_masking_by_bl(intree))
						
			treedata = parsed
			handle = StringIO(treedata)
			t = Phylo.read(handle, "newick")
		
			tips = []
			for leaf in t.get_terminals(): tips.append(leaf.name)
				
			if len(tips) != len(set(tips)): 
				print "Clades not monophyletic"
			
			else:
				
				with open(outDIR+base_name_tree[0]+".mt","w") as outfile:
					outfile.write(str(parsed)+";\n")
	
		
if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python modified_monophyly.py tree_DIR tree_file_ending clade_assignments outDIR"
		sys.exit(0)

	tree_DIR, tree_file_ending, clade_assignments, outDIR = sys.argv[1:]
	rename_tips(tree_DIR, tree_file_ending, clade_assignments, outDIR)



