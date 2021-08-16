"""
Taxon duplication? --No--> output one-to-one orthologs
		|
	   Yes
	    |
Outgroup present? --No--> ignore this homolog
		|
	   Yes
	    |
Outgroup taxon duplication? --Yes--> ignore this homolog
		|
	    No
	    |
Outgroup monophyletic? --No--> ignore this homolog
		|
	   Yes
	    |
Infer orthologs by using monophyletic, non-repeating outgroups
       
If not to output 1-to-1 orthologs, for example, already analysed these
set OUTPUT_1to1_ORTHOLOGS to False

"""

import phylo3,newick3,os,sys

OUTPUT_1to1_ORTHOLOGS = True

"""
INGROUPS = ["Alchemilla_alpina","Alchemilla_colura","Alchemilla_cryptantha","Alchemilla_decumbens","Alchemilla_ellenbeckii","Alchemilla_elongata","Alchemilla_fischeri","Alchemilla_fissa","Alchemilla_flabellata","Alchemilla_glabra","Alchemilla_haumanii","Alchemilla_heptagona","Alchemilla_hildebrandtii","Alchemilla_indivisa","Alchemilla_japonica","Alchemilla_johnstonii","Alchemilla_kiwuensis","Alchemilla_lapeyrousii","Alchemilla_microbetula","Alchemilla_pentaphyllea","Alchemilla_plicata","Alchemilla_roccatii","Alchemilla_rutenbergii","Alchemilla_schizophylla","Alchemilla_splendens","Alchemilla_stuhlmannii","Alchemilla_subnivalis","Alchemilla_subsericea","Alchemilla_triphylla","Alchemilla_woodii","Alchemilla_xanthochlora","Aphanes_arvensis","Aphanes_occidentalis","Aphanes_sp","Aphanes_cotopaxiensis","Aphanes_microcarpa","Aphanes_neglecta","Lachemilla_angustata","Lachemilla_equisetiformis","Lachemilla_llanganatensis","Lachemilla_moritziana","Lachemilla_pringlei","Lachemilla_ramosissima","Lachemilla_sibbaldiifolia","Lachemilla_standleyi","Lachemilla_velutina","Alchemilla_argyrophylla","Alchemilla_mollis","Aphanes_australis","Lachemilla_andina","Lachemilla_aphanoides","Lachemilla_barbata","Lachemilla_diplophylla","Lachemilla_erodiifolia","Lachemilla_galioides","Lachemilla_hirta","Lachemilla_hispidula","Lachemilla_jamesonii","Lachemilla_mandoniana","Lachemilla_nivalis","Lachemilla_orbiculata","Lachemilla_pectinata","Lachemilla_pinnata","Lachemilla_polylepis","Lachemilla_procumbens","Lachemilla_tanacetifolia","Lachemilla_verticillata","Lachemilla_vulcanica","Chamaerhodos_erecta","Comarum_palustre","Dasiphora_fruticosa","Drymocallis_arguta","Farinopsis_salesoviana","Sibbaldia_procumbens","Sibbaldianthe_adpressa","Sibbaldianthe_bifurca","Potaninia_mongolica","Fragaria_vesca","Drymocallis_glandulosa"]
OUTGROUPS = ["Rosa_woodsii","Potentilla_indica","Sanguisorba_menziesii"]
"""
INGROUPS = ["A","B","C","E","F"]
OUTGROUPS = ["D"]
"""

INGROUPS = ["Diabelia_serrata_E123","Diabelia_sanguinea_E301","Diabelia_ionostachya_stenophylla_E306","Diabelia_ionostachya_wenzhouensis_E204","Diabelia_spathulata_spathulata_E198","Abelia_macrotera_E110","Abelia_uniflora_E51","Abelia_forrestii_E37","Abelia_chinensis_ionandra_E30","Abelia_chinensis_E206","Kolkwitzia_amabilis_E9","Vesalea_occidentalis_E96","Vesalea_mexicana_E93","Vesalea_coriacea_E89","Vesalea_coriacea_E284","Vesalea_floribunda_E92","Linnaea_borealis_E14","Linnaea_borealis_E59","Dipelta_floribunda_E56","Dipelta_floribunda_E57","Dipelta_floribunda_E55","Scabiosa_canescens_E223","Scabiosa_tschiliensis_E21","Dipsacus_japonicus_E23","Valeriana_officinalis_E27","Valeriana_urticifolia_scorpioides_E219","Centranthus_ruber_E220","Valerianella_dentata_E217","Zabelia_biflora_E100","Zabelia_integrifolia_E15","Zabelia_dielsii_E108","Zabelia_dielsii_E286","Zabelia_triflora_E276","Morina_longifolia_E20","Acanthocalyx_alba_E19","Lonicera_korolkowii_E212","Lonicera_ligustrina_pileata_E74","Lonicera_confusa_E193","Lonicera_arizonica_E269","Symphoricarpos_orbiculatus_E237","Heptacodium_miconioides_E28","Weigela_florida_E99","Diervilla_lonicera_E331"]
OUTGROUPS = ["Viburnum_opulus_americanum_E162","Sambucus_williamsii_E208","Sambucus_nigra_E207"]
"""
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

	
def prune_paralogs_from_rerooted_homotree(root):
	if len(get_front_names(root)) == len(set(get_front_names(root))):
		return root #no pruning needed
	#check for duplications at the root first
	#one or two of the trifurcating root clades are ingroup clades
	node0,node1,node2 = root.children[0],root.children[1],root.children[2]
	out0,out1,out2=len(get_front_outgroup_names(node0)),len(get_front_outgroup_names(node1)),len(get_front_outgroup_names(node2))
	if out0 == 0 and out1 == 0: #0 and 1 are the ingroup clades
		name_set0 = set(get_front_names(node0))
		name_set1 = set(get_front_names(node1))
		if len(name_set0.intersection(name_set1)) > 0:
			if len(name_set0) > len(name_set1): #cut the side with less taxa
				root.remove_child(node1)
				node1.prune()		
			else:
				root.remove_child(node0)
				node0.prune()
	elif out1 == 0 and out2 == 0: #1 and 2 are the ingroup clades
		name_set1 = set(get_front_names(node1))
		name_set2 = set(get_front_names(node2))
		if len(name_set1.intersection(name_set2)) > 0:
			if len(name_set1) > len(name_set2): #cut the side with less taxa
				root.remove_child(node2)
				node2.prune()		
			else:
				root.remove_child(node1)
				node1.prune()
	elif out0 == 0 and out2 == 0: #0 and 2 are the ingroup clades
		name_set0 = set(get_front_names(node0))
		name_set2 = set(get_front_names(node2))
		if len(name_set0.intersection(name_set2)) > 0:
			if len(name_set0) > len(name_set2): #cut the side with less taxa
				root.remove_child(node2)
				node2.prune()		
			else:
				root.remove_child(node0)
				node0.prune()
	while len(get_front_names(root)) > len(set(get_front_names(root))):
		for node in root.iternodes(order=0): #PREORDER, root to tip
			if node.istip or node == root: continue
			child0,child1 = node.children[0],node.children[1]
			name_set0 = set(get_front_names(child0))
			name_set1 = set(get_front_names(child1))
			if len(name_set0.intersection(name_set1)) > 0:
				if len(name_set0) > len(name_set1): #cut the side with less taxa
					node.remove_child(child1)
					child1.prune()
				else:
					node.remove_child(child0)
					child0.prune()
				node,root = remove_kink(node,root) #no rerooting here
				break
	return root

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python python prune_paralogs_MO.py homoTreeDIR tree_file_ending minimal_taxa outDIR"
		sys.exit(0)
	
	inDIR = sys.argv[1]+"/"
	tree_file_ending = sys.argv[2]
	MIN_TAXA = int(sys.argv[3])
	outDIR = sys.argv[4]+"/"

	for i in os.listdir(inDIR):
		if not i.endswith(tree_file_ending): continue
		print i
		
		#read in the tree and check number of taxa
		outID = outDIR+get_clusterID(i)
		with open(inDIR+i,"r") as infile:
			 intree = newick3.parse(infile.readline())
		curroot = intree
		names = get_front_names(curroot)
		num_tips,num_taxa = len(names),len(set(names))
		if num_taxa < MIN_TAXA:
			continue #not enough taxa
		
		#If the homolog has no taxon duplication, no cutting is needed
		if num_tips == num_taxa:
			if OUTPUT_1to1_ORTHOLOGS:
				os.system("cp "+inDIR+i+" "+outID+".1to1ortho.tre")
		else:	
			#now need to deal with taxon duplications
			#check to make sure that the ingroup and outgroup names were set correctly
			for name in names:
				if name not in INGROUPS and name not in OUTGROUPS:
					print "check name",name
					sys.exit()
			outgroup_names = get_front_outgroup_names(curroot)
			
			#if no outgroup at all, do not attempt to resolve gene duplication
			if len(outgroup_names) == 0:
				print "duplicated taxa in unrooted tree"
				
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
					ortho = prune_paralogs_from_rerooted_homotree(curroot)
					if len(set(get_front_names(curroot))) >= MIN_TAXA:
						with open(outID+".ortho.tre","w") as outfile:
							outfile.write(newick3.tostring(ortho)+";\n")
					else: print "not enough taxa after pruning"
				else: print "outgroup non-monophyletic"
			
