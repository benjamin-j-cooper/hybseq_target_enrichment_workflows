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

OUTPUT_1to1_ORTHOLOGS = False

INGROUPS = ["Spergularia_heldreichii_PAFTOL_4935",
"Petrorhagia_dubia_PAFTOL_4936",
"Dianthus_fruticosus_PAFTOL_6029",
"Paronychia_argentea_PAFTOL_19651",
"Paronychia_kapela_PAFTOL_19653",
"Herniaria_glabra_PAFTOL_19655",
"Drymaria_cordata_PAFTOL_19657",
"Polycarpaea_repens_PAFTOL_19659",
"Loeflingia_hispnica_PAFTOL_19661",
"Pycnophyllum_bryoides_PAFTOL_19663",
"Cardionema_ramosissima_PAFTOL_19665",
"Gypsophila_repens_PAFTOL_19667",
"Acanthophyllum_mucronatum_PAFTOL_19669",
"Spergularia_heldreichii_PAFTOL_4935",
"Petrorhagia_dubia_PAFTOL_4936",
"Dianthus_fruticosus_PAFTOL_6029",
"Paronychia_argentea_PAFTOL_19651",
"Paronychia_kapela_PAFTOL_19653",
"Herniaria_glabra_PAFTOL_19655",
"Drymaria_cordata_PAFTOL_19657",
"Polycarpaea_repens_PAFTOL_19659",
"Loeflingia_hispnica_PAFTOL_19661",
"Pycnophyllum_bryoides_PAFTOL_19663",
"Cardionema_ramosissima_PAFTOL_19665",
"Gypsophila_repens_PAFTOL_19667",
"Acanthophyllum_mucronatum_PAFTOL_19669",
"Spergularia_heldreichii_PAFTOL_4935",
"Petrorhagia_dubia_PAFTOL_4936",
"Dianthus_fruticosus_PAFTOL_6029",
"Paronychia_argentea_PAFTOL_19651",
"Paronychia_kapela_PAFTOL_19653",
"Herniaria_glabra_PAFTOL_19655",
"Drymaria_cordata_PAFTOL_19657",
"Polycarpaea_repens_PAFTOL_19659",
"Loeflingia_hispnica_PAFTOL_19661",
"Pycnophyllum_bryoides_PAFTOL_19663",
"Cardionema_ramosissima_PAFTOL_19665",
"Gypsophila_repens_PAFTOL_19667",
"Acanthophyllum_mucronatum_PAFTOL_19669",
"Moehringia_muscosa_PAFTOL_19671",
"Pseudostellaria_heterophylla_PAFTOL_19673",
"Stellaria_graminea_PAFTOL_19675",
"Mesostemma_kotschyanum_subsp_afghanicum_PAFTOL_19677",
"Cerastium_arvense_PAFTOL_19679",
"Sagina_procumbens_PAFTOL_19681",
"Minuartia_recurva_PAFTOL_19683",
"Bufonia_macropetala_PAFTOL_19685",
"Drypis_spnosa_PAFTOL_19687",
"Habrosia_spnuliflora_PAFTOL_19689",
"Scleranthus_annuus_PAFTOL_19691",
"Mcneillia_graminifolia_subsp_rosanoi_PAFTOL_19693",
"Arenaria_serpyllifolia_PAFTOL_19695",
"Stellaria_Holostea_PAFTOL_19697",
"Odontostemma_leucasterium_PAFTOL_19699",
"Arenaria_globiflora_PAFTOL_19701",
"Sabulina_tenuifolia_PAFTOL_19703",
"Sabulina_decandra_PAFTOL_19705",
"Stellaria_filiformis_GAP_22101",
"Colobanthus_apetalus_GAP_22145",
"Sagina_maritima_GAP_22425",
"Scleranthus_pungens_GAP_22461",
"Colobanthus_apetalus_GAP_23267",
"Sagina_maritima_GAP_23321",
"Scleranthus_pungens_GAP_23323",
"Stellaria_filiformis_GAP_23325",
"Lepyrodiclis_holosteoides_PAFTOL_25171",
"Polytepalum_angolense_PAFTOL_25173",
"Augustea_suffruticosa_PAFTOL_25177",
"Cerastium_lanceolatum_PAFTOL_25201",
"Cometes_abyssinica_PAFTOL_25297",
"Graecobolanthus_graecus_PAFTOL_25313",
"Thylacosprmum_caesptosum_PAFTOL_25331",
"Rhodalsine_geniculata_PAFTOL_25341",
"Bolanthus_filicaulis_PAFTOL_25343",
"Shivparvatia_glanduligera_PAFTOL_25397",
"Gymnocarpos_decander_PAFTOL_25415",
"Silene_alpestris_PAFTOL_25435",
"Stellaria_monosprma_PAFTOL_25461",
"Achyronychia_cooperi_PAFTOL_25481",
"Stellaria_obtusa_PAFTOL_25483",
"Stipulicida_setacea_PAFTOL_25485",
"Arenaria_patula_PAFTOL_25487",
"Stellaria_filiformis_GAP_26429",
"Colobanthus_apetalus_GAP_26443",
"Sagina_maritima_GAP_26633",
"Scleranthus_pungens_GAP_26669"]
OUTGROUPS = ["beta_vulgaris_CDS_outgroup",
"mesembryanthemum_crystallinum_CDS_outgroup",
"simmondsia_chinensis_CDS_outgroup",
"spinacia_oleracea_CDS_outgroup"]


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
		print ("python python prune_paralogs_MO.py homoTreeDIR tree_file_ending minimal_taxa outDIR")
		sys.exit(0)

	inDIR = sys.argv[1]+"/"
	tree_file_ending = sys.argv[2]
	MIN_TAXA = int(sys.argv[3])
	outDIR = sys.argv[4]+"/"

	for i in os.listdir(inDIR):
		if not i.endswith(tree_file_ending): continue
		print (i)

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
					print ("check name",name)
					sys.exit()
			outgroup_names = get_front_outgroup_names(curroot)

			#if no outgroup at all, do not attempt to resolve gene duplication
			if len(outgroup_names) == 0:
				print ("duplicated taxa in unrooted tree")

			#skip the homolog if there are duplicated outgroup taxa
			elif len(outgroup_names) > len(set(outgroup_names)):
				print ("outgroup contains taxon repeats")

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
					else: print ("not enough taxa after pruning")
				else: print ("outgroup non-monophyletic")
