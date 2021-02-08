"""
Script to root trees with multiple taxa. Modify from pruning MO. It will check that there is no duplicated taxa in the outgroup and require 'monophyletic' outgroups (i.e outgroups nested in the  ingroup taxa). Mainly to use with Phyparts. 
Outgroup names need to be provided in text file with each taxa on a line
"""
import phylo3,newick3,os,sys

OUTGROUPS = []

def check_outgroup_file(outgroup_list):
	for line in open(outgroup_list,'r').readlines():
		OUTGROUPS.append(line.strip())
	print ("outgroups:")
	print (OUTGROUPS)

#given tip label, return taxon name identifier


def get_name(label):
	return label.split("@")[0]

def get_clusterID(filename):
	return filename.split(".")[0]

def get_front_labels(node):
	leaves = node.leaves()
	return [i.label for i in leaves]

def get_back_labels(node,root):
	all_labels = get_front_labels(root)
	front_labels = get_front_labels(node)
	return set(all_labels) - set(front_labels) #labels do not repeat
	
def get_front_names(node): #may include duplicates
	labels = get_front_labels(node)
	return [get_name(i) for i in labels]

def get_front_outgroup_names(node):
	names = get_front_names(node)
	return [i for i in names if i in OUTGROUPS]

def get_back_names(node,root): #may include duplicates
	back_labels = get_back_labels(node,root)
	return [get_name(i) for i in back_labels]

def remove_kink(node,curroot):
	if node == curroot and curroot.nchildren == 2:
		#move the root away to an adjacent none-tip
		if curroot.children[0].istip: #the other child is not tip
			curroot = phylo3.reroot(curroot,curroot.children[1])
		else: curroot = phylo3.reroot(curroot,curroot.children[0])
	#---node---< all nodes should have one child only now
	length = node.length + (node.children[0]).length
	par = node.parent
	kink = node
	node = node.children[0]
	#parent--kink---node<
	par.remove_child(kink)
	par.add_child(node)
	node.length = length
	return node,curroot
	
#check if outgroups are monophyletic and non-repeating and reroot
#otherwise return None
def reroot_with_monophyletic_outgroups(root):
	lvs = root.leaves()
	outgroup_matches = {} #key is label, value is the tip node object
	#Since no taxon repeat in outgroups name and leaf is one-to-one
	outgroup_labels = []
	for leaf in lvs:
		label = leaf.label
		name = get_name(label)
		if name in OUTGROUPS:
			outgroup_matches[label] = leaf
			outgroup_labels.append(leaf.label)
	if len(outgroup_labels) == 1: #one single outgroup
		#cannot reroot on a tip so have to go one more node into the ingroup
		new_root = outgroup_matches[outgroup_labels[0]].parent
		return phylo3.reroot(root,new_root)
	else: #has multiple outgroups. Check monophyly and reroot
		newroot = None
		for node in root.iternodes():
			if node == root: continue #skip the root
			front_names = get_front_names(node)
			back_names = get_back_names(node,root)
			front_in_names,front_out_names,back_in_names,back_out_names = 0,0,0,0
			for i in front_names:
				if i in OUTGROUPS: front_out_names += 1
				else: front_in_names += 1
			for j in back_names:
				if j in OUTGROUPS: back_out_names += 1
				else: back_in_names += 1
			if front_in_names==0 and front_out_names>0 and back_in_names>0 and back_out_names==0:
				newroot = node #ingroup at back, outgroup in front
				break
			if front_in_names>0 and front_out_names==0 and back_in_names==0 and back_out_names>0:
				newroot = node.parent #ingroup in front, outgroup at back
				break
		if newroot != None:
			return phylo3.reroot(root,newroot)
		else: return None

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python root_trees_mulitple_outgroup_MO.py inDIR tree_file_ending outDIR outgroup_list"
		sys.exit(0)
	
	inDIR = sys.argv[1]+"/"
	tree_file_ending = sys.argv[2]
	outDIR = sys.argv[3]+"/"
	outgroup_list = sys.argv[4]

	check_outgroup_file(outgroup_list)

	for i in os.listdir(inDIR):
		if not i.endswith(tree_file_ending): continue
		print i
		
		outID = outDIR+get_clusterID(i)
		with open(inDIR+i,"r") as infile:
			 intree = newick3.parse(infile.readline())
		curroot = intree
		
		outgroup_names = get_front_outgroup_names(curroot)
			
		#if no outgroup at all, do not attempt to resolve gene duplication
		if len(outgroup_names) == 0:
			print "no outgroups present"
				
		#skip the homolog if there are duplicated outgroup taxa
		elif len(outgroup_names) > len(set(outgroup_names)): 
			print "outgroup contains taxon repeats"
				
		else: #at least one outgroup present and there's no outgroup duplication
			if curroot.nchildren == 2: #need to reroot
				temp,curroot = remove_kink(curroot,curroot)
			curroot = reroot_with_monophyletic_outgroups(curroot)
			#only return one tree after prunning
			if curroot != None:
				with open(outID+".reroot","w") as outfile:
					outfile.write(newick3.tostring(curroot)+";\n")
			else: print "outgroup non-monophyletic"
			
