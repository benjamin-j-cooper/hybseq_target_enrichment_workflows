# Target enrichment orthology


#### This repository contains scripts for orthology in inference from target enrichment data (e.g. Hyb-Seq) and relies mainly in scripts from [Phylogenomic dataset construction respository](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/master/) plus some new ones.

#### If using the scrips from this repository you must cite

Morales-Briones, D.F., B. Gehrke, H. Chien-Hsun Huang, A. Liston, M. Hong. H.E. Marx, D.C. Tank & Y. Yang. Analysis of paralogs in target enrichment data pinpoints multiple ancient polyploidy events in Alchemilla s.l. (Rosaceae). [bioRxiv 2020.08.21.261925; doi: https://doi.org/10.1101/2020.08.21.261925 ](https://doi.org/10.1101/2020.08.21.261925)

Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution. [doi:10.1093/molbev/msu245](https://doi.org/10.1093/molbev/msu245)


### Dependencies needed to run the scripts. 

[TreeShrink](https://github.com/uym2/TreeShrink) It works now with Version 1.3.2 (older versions won't work)

[RAxML](https://github.com/stamatak/standard-RAxML) Version 8.2.11  (newer versions should work)

[Phyx](https://github.com/FePhyFoFum/phyx)

[MAFFT](https://mafft.cbrc.jp/alignment/software/) Version 7.307 (newer versions should work)

[MACSE](https://bioweb.supagro.inra.fr/macse/index.php?menu=releases) Version 2.0.3 (newer versions should work)


## Step 1: Format paralog fasta files

#### **This example is for the output of 'paralog_investigator' of [HybPiper](https://github.com/mossmatters/HybPiper/wiki/Paralogs)**

To use the script from [Phylogenomic dataset construction respository](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/master/) an "@" symbol needs to be added in the fasta headers to identify paralogs of the same sample.

The fasta files *.paralogs.fasta after running 'paralog_investigator' have a format from [SPAdes](http://cab.spbu.ru/software/spades/) like

\>Rosa_woodsii  
AGTC...  
\>Lachemilla_pinnata.0 NODE_2_length_1787_cov_230.693567,Fragaria-gene15996_1557_01,4,519,78.47,(+),1,1537  
ACGT.....  
\>Alchemilla_colura.main NODE_3_length_1706_cov_62.896426,Fragaria-gene15996_1557_01,0,517,81.43,(-),1552,1  
ACCC....  
\>Alchemilla_colura.1 NODE_1_length_2101_cov_47.514174,Fragaria-gene15996_1557_01,0,519,79.11,(+),136,1687  
ACCG....  

##### To format the sequences run the next loop in the directory where fasta files are locaded: (Note: This command will overwrite the fasta files)

	for i in $(ls *.fasta); do
	sed -i -E 's/(>.+)\.(.+)\s(.+)/\1@paralog_\2/‘ $i
	sed -i '/>/ {/@/! s/$/@unique/}’ $i
	done 
	
The output of this will be

\>Rosa_woodsii@unique  
AGTC...  
\>Lachemilla_pinnata@paralog_0  
ACGT.....  
\>Alchemilla_colura@paralog_main  
ACCC....  
\>Alchemilla_colura@paralog_1  
ACCG....  

##### For details of SPAdes contigs see [HybPiper's paralogs description](https://github.com/mossmatters/HybPiper/wiki/Paralogs)


If other additional sequences are added to the fasta files (e.g. from reference genomes) make sure that those also have the "@" format.

I added reference sequences of *Fragaria vesca* as Fragaria_vesca@genome


## Step 2: Build homolog trees

I used Macse for DNA alignment. Macse does not have multithread running options, so I wrote individual bash files to run them in parallel.

##### To write bash files

	for filename in $(ls *.fasta)
	do
	echo java -jar ~/Apps/macse_v2.03.jar -prog alignSequences -seq $filename -out_NT $(cut -d'.' -f1 <<<"$filename").NT.aln -out_AA $(cut -d'.' -f1 <<<"$filename").AA.aln > $(cut -d'.' -f1 <<<"$filename").sh
	done  


##### To run alignment runs in parallel

	parallel -j 32 bash ::: *.sh

If there are frame shifts in the alignments, Macse will replace the shifted codons with "!" and will cause problems with RAxML or IQtree. 


##### To replace "!" codon with gaps

	python remove_shifted_codons_from_macse.py <alignment directory> <fasta file extension> <output directory> <nt or aa>
		
	
##### Alternatively, Mafft can be used for the alignment

 	python mafft_wrapper.py <fasta files directory> <fasta file extension> <# of threads> <dna or aa>
 	
 	
##### Trim alignments with Phyx

	python pxclsq_wrapper.py <alignment directory > <mininum column occupancy> <dna or aa>
	
The output will files with extension .aln-cln


##### Infer ML trees with RAxML. This will infer trees with GTRGAMMA and 100 bootstrap replicates (different model and # of bs replicates can be modify in the script).  

	python raxml_bs_wrapper.py <aln-cln files directory> <# of threads> <dna or aa>
	
	
##### If no bs replicates are needed use. 
	
	python raxml_wrapper.py <aln-cln files directory> <# of threads> <dna or aa>
	
	
##### Mask both mono- and (optional) paraphyletic tips that belong to the same taxon. Keep the tip that has the most un-ambiguous charactors in the trimmed alignment.

	python mask_tips_by_taxonID_transcripts.py <tre files directory> <aln-cln files directory> mask_paraphyletic <y or n>
	

##### Alternatively, to remove only monophyletic tips and keep the sequence with the shortest terminal branch length.

	python mask_tips_by_taxonID_genomes.py <tre files directory>
	

##### Trim spurious tips with [TreeShrink](https://github.com/uym2/TreeShrink)

	python tree_shrink_wrapper.py <input directory> <tree file extension> <quantile>

It outputs the tips that were trimmed in the file .txt and the trimmed trees in the files .tt. You would need to test different quantiles to see which one fit better your data. 


##### These are you final homologs trees. If you want to repeat the tree inference, masking, and trimming or just re-infer homologs trees without removed tips. Write fasta files from trimmed trees. 

	python write_homolog_fasta_from_multiple_aln.py <original "@" formated fasta files directory> <trimmed trees directory> <fasta file extension> <tree file extension> <output directory>


## Step 3: Paralogy pruning to infer orthologs. Use one of the following


#### For details of orthology inference refer to Yang and Smith (2014) [doi:10.1093/molbev/msu245](https://doi.org/10.1093/molbev/msu245)


##### 1to1: only look at homologs that are strictly one-to-one. No cutting is carried out.
	
	python filter_1to1_orthologs.py <final homologs directory> <tree file extension> <minimal number of taxa> <output directory>


##### MI: prune by maximum inclusion.  Set OUTPUT_1to1_ORTHOLOGS to False if wish only to ouput orthologs that is not 1-to-1, for example, when 1-to-1 orthologs have already been analyzed in previous steps.

	python prune_paralogs_MI.py <final homologs directory> <tree file extension> <relative long tip cutoff> <absolute long tip cutoff> <minimal number of taxa> <output directory>


##### MO: prune by using homologs with monophyletic, non-repeating outgroups, reroot and cut paralog from root to tip. If no outgroup, only use those that do not have duplicated taxa. Change the list of ingroup and outgroup names first. Set OUTPUT_1to1_ORTHOLOGS to False if wish only to ouput orthologs that is not 1-to-1

	python prune_paralogs_MO.py <final homologs directory> <tree file extension> <minimal number of taxa> <output directory>


##### RT: prune by extracting ingroup clades and then cut paralogs from root to tip. If no outgroup, only use those that do not have duplicated taxa. Compile a list of ingroup and outgroup taxonID, with each line begin with either "IN" or "OUT", followed by a tab, and then the taxonID.

	python prune_paralogs_RT.py <final homologs directory> <tree file extension> <output directory>  <minimal number of taxa> <ingroup and outgroup taxonIDs table>


##### Or alternatively, if the input homolog tree is already rooted

	python prune_paralogs_from_rooted_trees.py <final homologs directory> <tree file extension> <minimal number of taxa> <output directory>


## Step 4: Visualize matrix occupancy stats, write final fasta files and build supermatrix (optional)

	python ortholog_occupancy_stats.py <orthologs directory>


##### Read in and rank number of taxa per ortholog from highest to lowest. Plot in R the ranked number of taxa per ortholog

	a <- as.numeric(read.table("ortho_stats")[1])
	a <- sort(a, decreasing=TRUE)
	pdf(file="taxon_occupancy.pdf")
	plot(a, type="l", lwd=3, ylab="Number of Taxa in Each Ortholog")
	dev.off()


##### Check taxon_stats to see if any taxa have unusally low number of genes in the orthologs and decide the minimum number of taxa for the supermatrix. 

##### Write new and final fasta files from ortholog trees

	python write_ortholog_fasta_from_multiple_aln.py <original "@" formated fasta files directory> <orthologs directory> <fasta file extension> <orthologs tree file extension> <output directory>
	
##### With those you can align, clean, infer trees from final orthologs, and infer species trees using your preferred tools

 
##### If a supermatrix of clean alignment is needed. Once you have clean alignments, choose the minimal cleaned alignment length and minimal number of taxa of each ortholog to include in the supermatrix

	python concatenate_matrices_phyx.py <aln-cln files directory> <minimum number of sites> <minimum number of taxa> <output directory>

This will output a list of cleaned orthology alignments that passed the filter, a summary of taxon matrix occupancies to check whether any taxon is under represented, and a concatenated matrix in fasta, phylip, and nexus format as well as a partition file for RAxML or IQtree.






