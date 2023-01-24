# Target enrichment orthology


#### This repository contains instructions and scripts for running the hybpiper program with orthology inference.
#### It is intended to be used with target enrichment data (e.g. Hyb-Seq), however there are options to include data from genomes (see [plastome wiki page](https://github.com/benjamin-j-cooper/hybseq_target_enrichment_workflows/wiki/Target-enrichment-plastome-analysis) ).

#### These instructions rely mainly on scripts from [Phylogenomic dataset construction respository](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/master/) and [target_enrichment_orthology](https://bitbucket.org/dfmoralesb/target_enrichment_orthology/src/master/) plus some new ones.

#### If using the scrips from this repository you must cite

Morales-Briones, D.F., B. Gehrke, H. Chien-Hsun Huang, A. Liston, M. Hong. H.E. Marx, D.C. Tank & Y. Yang.  Analysis of paralogs in target enrichment data pinpoints multiple ancient polyploidy events in *Alchemilla* s.l. (Rosaceae). [ Systematic Biology, syab032](https://doi.org/10.1093/sysbio/syab032)

Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution. [doi:10.1093/molbev/msu245](https://doi.org/10.1093/molbev/msu245)


# Dependencies

python2.7 and python 3.9 or later 
Because both python2 and python3 are currently required for this workflow, I highly recommend using Conda to manage you python installations. this makes it easy to switch between python2 and python3 environments. 

[GNUparallel](https://www.gnu.org/software/parallel/) 2018 version or later

[hybpiper](https://github.com/mossmatters/HybPiper) Legacy version 1.3 (To run this pipeline with Hybpiper 2, see the [target enrichment orthology with hybpiper 2 wiki](https://github.com/benjamin-j-cooper/hybseq_target_enrichment_workflows/wiki/Target-Enrichment-Orthology-with-Hybpiper-2))

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.5 or later

[multiqc](https://multiqc.info) v1.12 or later

[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.39 or later

[samtools](https://www.htslib.org) v1.15 or later

[TreeShrink](https://github.com/uym2/TreeShrink) 1.3.2 or later

[RAxML](https://github.com/stamatak/standard-RAxML) 8.2.11 or later

[Phyx](https://github.com/FePhyFoFum/phyx) v0.1 or later

[MAFFT](https://mafft.cbrc.jp/alignment/software/) 7.307 or later

[OMM_MACSE](https://bioweb.supagro.inra.fr/macse/index.php?menu=releases) v10.02 or later

### To run hybpiper, you will need a target file. I recommend making a custom target file with samples from the family of your study. 
See the [hybpiper help page](https://github.com/mossmatters/HybPiper/wiki/Troubleshooting,-common-issues,-and-recommendations) for other suggestions.		

# Part 1: data preparation
##  1.1. Download raw reads
### 1.1.1. In your working directory, make directories:		

	mkdir /{project_directory}/raw_reads
	mkdir -p /{project_directory}/raw_reads/fastqc
	mkdir -p /{project_directory}/trimmed_reads/unpaired_trimmed
	mkdir -p /{project_directory}/trimmed_reads/trimm_stats
	mkdir -p /{project_directory}/trimmed_reads/fastqc
	mkdir -p /{project_directory}/hybpiper_output/orthology
	mkdir -p /{project_directory}/hybpiper_output/outgroup_genomes
	

### 1.1.2. Make sample files
first cd in into the raw_reads directory

	cd raw_reads
		
I wrote a custom script to do this next part for the PAFTOL dataset [make_sample_files_from_drive.py](https://github.com/benjamin-j-cooper/phylo_genomic_utilities), it produces text files needed for downloading, renaming, and calling samples
You can edit this script to make files for your project if it is not part of the PAFTOL project.		
To run the script:

	python make_sample_files_from_drive.py -f {family}

#### If you used the above script, skip to step 1.1.4
#### 1.1.3. if are not using the script, you will need to manually make three files:
1.  download_string.txt
	each line of file is wget command to download a read file
	save this in /{project_directory}/raw_reads/
	
	

2.  filename_key.txt 
	This is a tab delimited file with first column containing original file name and second column containing new file name (genus_species_subspecies_sampleID.fq.gz)		
	/{project_directory}/raw_reads/
	
	
	
3.  sample_names.txt
	This is a list with the sample names you want for each sample as a line, I chose the root name from the new file names
	genus_species_subspecies_sampleID	
	save this in /{project_directory}/hybpiper_output/
	
	
	
### 1.1.4. Download raw reads into the raw_reads folder from sftp link or other source using the download_string.txt created in earlier and batch script download.sh
open screen and run script from raw_reads directory:
	
	screen -S download
	parallel wget {} :::: download_string.txt


#### 1.1.5. Once you download you will need to rename files:
I parallelized the mv function in linux to rename downloaded files using the filename_key.txt: 

	parallel --colsep '\t' mv {1} {2} :::: filename_key.txt

## 1.2 QC raw data
### 1.2.1 QC raw data with fastqc	
from raw_reads directory, make a txt file with names of all the .fq files to QC, then run fastqc

	ls *fastq.gz > ./fastqc/fastq_files.txt		
	
then run fastqc		

	while read i; do fastqc -f fastq -t 6 $i -o /{project_directory}/raw_reads/fastqc/ --noextract; done < ./fastqc/fastq_files.txt

#### 1.2.2. then from fastqc directory run multiqc
multiqc compiles fastqc reports into one larger report for easy veiwing
	
	cd fastqc
	multiqc . -n multiqc_raw_reads_report.html

### 1.3 Remove adapters and quality filter
Here I use an Illumina adapter sequence file called "TruSeq_adapters.fa" it is publicly available [here]() or you can use your own tailored to your project.

#### 1.3.1 In a text editor, make trimmomatic_GEN.sh script
save it in the raw_reads directory, it will generate a script file for every read file (adjust trimmomatic settings as desired)
Here is the code that I used for the script:

	for f1 in *1.fastq.gz
	do
	f2=${f1%%1.fastq.gz}"2.fastq.gz"
	echo "java -jar -Xmx10g  /{apps_directory}/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 $f1 $f2 ../${f1%%.fastq.gz}_paired.fq.gz ../unpaired_trimmed/${f1%%.fastq.gz}_unpaired.fq.gz ../${f2%%.fastq.gz}_paired.fq.gz ../unpaired_trimmed/${f2%%.fastq.gz}_unpaired.fq.gz ILLUMINACLIP:/target_enrichment_orthology/files/TruSeq_adapters.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:25 2> ../trimm_stats/trimmomatic_stats_${f1%%1.fastq.gz}.txt" > "$f1.trimm.sh"
	done

#### 1.3.2 activate script		

	chmod +x trimmomatic_GEN.sh

#### 1.3.3 call trimmomatic_GEN.sh to make trimmomatic scripts for each sample

	./trimmomatic_GEN.sh

#### 1.3.4 activate all these scripts

	chmod +x *.trim.sh

#### 1.3.5 Run all individual trimmomatic scripts using GNUparallel 
note that the trimmomatic scripts will put the paired and trimmed files up one level in the folder hierarchy leaving the raw reads behind in the raw_reads folder and placing the other parts of the trimmomatic output in appropriate folders.

	parallel -j 12 bash {} ::: *.fastq.gz.trim.sh

#### 1.3.6 unzip trimmed and unpaired_trimmed files
	
	parallel gunzip {} ::: *.fq.gz
	
	cd unpaired_trimmed
	parallel gunzip {} ::: *.fq.gz

#### 1.3.7 QC trimmed reads with fastqc
from __trimmed_reads__ directory, make a txt file with names of all the .fq files to QC, then run fastqc

	ls *paired.fq > ./fastqc/fastq_files.txt
	while read i; do fastqc -f fastq -t 6 $i -o /{project_directory}/trimmed_reads/fastqc/ --noextract; done < ./fastqc/fastq_files.txt

#### 1.3.8 then run multiqc
cd into the fastqc directory and run multiqc

	cd fastqc
	multiqc . -n multiqc_trimmed_reads_report.html

# Part 2: hybpiper
### 2.1 Run hybpiper
	cd /{project_directory}/hybpiper_output/
	
To get paralogs with different percentages you need to first run a full HyPipper assembly
make hybpiper script and save as hybpiper.sh
notice that it is run in parrallel using the sample_name.txt file created earlier and specifying the number of cpu's with --cpu

	while read i
	do
	echo "python /{apps_directory}/hybPiper/reads_first.py --prefix $i --bwa --cpu 10 -b target_file.fasta -r /{project_directory}/trimmed_reads/$i*_paired.fq --unpaired /{project_directory}/trimmed_reads/unpaired_trimmed/${i%%}_cat_unpaired.fq" > "/{project_directory}/hybpiper_output/$i.hybpiper.sh"
	done < /{project_directory}/hybpiper_output/sample_names.txt

then activate the script

	chmod +x hybpiper.sh

and run

	./hybpiper.sh

### 2.2 Rerun hybpiper with modified version for less stringent paralog 
The percentage for calling paralogs is in line 459 of the 'exonerate_hits.py' script. The original line is "longhits = [x > 0.75*protlength for x in hitlengths]"
I changed the percentage in line 459 in the modified version of hybpiper in this repository to > 50%
Now run the modified version of hybPiper from this repository in the same hybpiper_output folder but skipping the BLAST, read distribution, mapping and assembly steps:

make second hybpiper script and save as hybpiper_50.sh

	while read i
	do
	echo "python /target_enrichment_orthology/scripts/reads_first_50.py --prefix $i --no-blast --no-distribute --no-assemble --bwa --cpu 10 -b target_file.fasta -r /{project_directory}/trimmed_reads/$i*_paired.fq --unpaired /{project_directory}/trimmed_reads/unpaired_trimmed/${i%%}_cat_unpaired.fq" > "/{project_directory}/hybpiper_output/$i.hybpiper_50.sh"
	done < /{project_directory}/hybpiper_output/sample_names.txt

activate script and run
	
	chmod +x hybpiper_50.sh
	./hybpiper.sh

Rerunning this step is pretty fast because it uses the previous run results. So, you need to run the full HybPiper assembly only once and only the output of exonerate will be modified. 

### 2.3 cleanup un-neaded intermediate files created by hybpiper

	while read i; do python /{apps_directory}/hybPiper/cleanup.py $i; done < sample_names.txt

### 2.4 Run stats on gene recovery and assembly:
nifty script to check enrichment efficiency

	while read i; do samtools flagstat $i/$i.bam; done < sample_names.txt

this is for making the gene_recovery_heatmap with the R script provided with the HybSeq wiki

	python /{apps_directory}/hybPiper/get_seq_lengths.py target_file.fasta sample_names.txt dna > gene_lengths.txt

save gene_lengths.txt to computer
see HybPiper documentation on GITHub for R script to generate heatmap

### 2.5 Make genelist.txt file from target_file.fasta
To do this, I made a copy of my target_file.fasta and opened the copy in a text editor on my local computer (I used BBedit, but you can use a text editor of your choosing).		
Using find and replace with RegEx, I broke each fasta header and sequence into capture groups. For my fasta headers: 

	>Arath-4471E01
	AACGTGGTTGAAGATAAAGAAAGGCTTGAGACTGCTAATACAGATTGGATGCATAAGTAC
	AAAGGATCTAGTAAGCTGATGCTCTTGCCCAAGAATACACAAGAG
	>NHUA-4471E01
	AATGTGATTGAAGATGAAGACGCGCTTCTCTTTGCAAACACTGATTGGATGCGAAAATAT
	AAGGGCTCCAGTAAGCTTTTGCTACAGCCTAGAACCACTCAAGAG
	>AJFN-4471E01
	AATGTTATACAGGATGAAGAGAAACTGAATACTGCAAACTCCGATTGGATGCGGAAATAC
	AAAGGCTCAAGTAAGCTTATGCTCCAACCTAGGAGCACCGAGGAG

My RegEx looked something like this, which is fairly generic and should work for most hybpiper formatted target files:

find:

	(>[A-Za-z]*-)([A-Za-z0-9]*\n)([A-Z\n]*)
	
replace:

	\2

Using Capture group \2 as the replace statement left me with a list of gene names in this format:

	4471E01
	4471E01
	4471E01

Because my target file had multiple samples per gene, I was then left with duplicate gene names. To easily remove duplicated I used the "process duplicate lines..." tool in the Text menu in BBedit. 
With duplicates removed from the the gene list, save this file as genelist.txt in your hybpiper_output directory

# Part 3
## 3.1 Make outgroups from genomes		
In your outgroup_genomes directory in your hybpiper_output directory, make a directory for each outgroup

	cd outgroup_genomes
	mkdir beta_vulgaris
	
Download the fasta files for the genome CDS's of your choosing into the outgroups' directory. In this example, I use [beta vulgaris](https://www.ncbi.nlm.nih.gov/nuccore/NC_059019.1?report=fasta) and downloaded the CDS manually onto my computer (you could use NCBI command line tools to do this directly to your remote workstation) and then I used filezilla to upload to my project directory.

Then, for each outgroup from the outgroups' named folder run [blast2fasta.py](https://github.com/benjamin-j-cooper/phylo_genomic_utilities) referencing the target file you used for hybpiper:

	python blast2fasta.py -i EL10_2_2.cdna.fa -o beta_vulgaris -q ../../caryophyllaceae_filtered_target_file.fasta -t nucl

Followed by [targets2hybpiper_directory](https://github.com/benjamin-j-cooper/phylo_genomic_utilities), referencing the genelist.txt file you created earlier in the hybpiper_output directory.

	python targets2hybpiper_directory.py -s beta_vulgaris_CDS_outgroup -g ../../genelist.txt		

*Make sure you at the chosen names of your outgroups' directory (-s) from this step to the sample_names.txt file located in the hybpiper_output directory.		
Finally, copy the newly created hybpiper-structured outgroup directory and subdirectories to your hybpiper_output root directory:

	cp -r beta_vulgaris_CDS_outgroup ../../		

Now you are ready to gather assemblies from hybpiper (including outgroups) by gene using hybpipers built in tools.

# Part 4 
### 4.1 paralog investigator
Now use 'paralog_investigator.py' and 'paralog_retriever.py' to get the alignments with the percentage that you want to test.
	
	while read i; do echo $i python /{apps_directory}/hybPiper/paralog_investigator.py  $i >> ./orthology/paralogs.txt; done < sample_names.txt

edit paralog_long.txt so there is only unique gene for each line. I did this in BB edit but you can do this in any text editor using find and replace:

	^\d+?\sparalogs\swritten\sfor\s 

and replace with nothing

then remove duplicates - In BBEdit, you can do this easily by going to "Text->Process Duplicated Lines"
then use find and replace in text editor to replace all \r with a single space to make a list of genes separated by spaces
save this file as paralogs_unique.txt in your /orthology/ directory

### 4.2 Paralog retriever
Run for select genes identified by paralog investigator (I exported the .txt file I made from the output screen of paralog investigator into R, removed the headers and then selected unique gene names with a custom script. the file is called "paralogs_50_curated" in my hybpiper directory):
Using a file with one gene per line (this is currently not working but would be an easy fix)

	parallel "python /{apps_directory}/hybPiper/paralog_retriever.py sample_names.txt {} > ./orthology/{}.paralogs.fasta" ::: paralogs_unique.txt 2> ./orthologs/paralog_table.txt

#using gene names (this is working):

	parallel "python /{apps_directory}/hybPiper/paralog_retriever.py sample_names.txt {} > ./orthology/{}.paralogs.fasta" ::: 6051 6295 4527 5664 5460 6439 6026 4802 6947 4724 6238 5821 6176 6641 6909 4757 6282 6538 6992 6447 5328 6494 5721 6303 6119 6000 6036 6780 4989 6462 5842 6373 5594 5942 6387 6779 6875 6068 5131 6528 7313 6527 5206 7628 5271 6483 5857 6274 7324 6460 5910 5355 6914 4954 6883 6732 5280 5449 5865 5434 5463 5162 6130 6550 6782 6003 5913 5357 5430 5980 6882 5802 6825 5273 5981 6034 7028 7371 6886 6533 5919 6412 5428 7333 6620 5941 5744 5343 6498 6946 5426 6791 6487 6532 6128 5536 5620 5188 5950 4951 5958 5469 6110 5333 6500 5933 7367 6559 6531 6148 7336 6164 6407 6299 6955 6401 6979 5660 6449 5220 4932 5578 6557 5177 7067 6954 6488 7325 6139 5422 7194 5840 5853 6492 6785 5347 5656 5116 5531 6284 7602 5940 6797 6705 6379 6969 6860 6962 7583 4889 6038 7029 6048 5528 5326 6913 5859 6227 7273 5960 6420 5257 7363 7024 6933 6198 6496 5464 5702 5299 4691 4992 6114 5945 5949 6016 6064 5398 5893 6738 6221 6631 5123 6649 6958 5562 5090 6978 5064 5894 5642 6506 7361 5034 7141 6459 5477 6733 6150 5974 6004 7111 6265 6393 5168 2> ./orthologs/paralog_table.txt

cd into orthology directory 

	cd orthology

remove empty genes files where no sequences were captured (this may resolve issues with RAxML later on)

	find ./*paralogs.fasta -maxdepth 1 -type f -empty -delete > empty_genes.txt

# Part 5 orthology

	cd orthology
	mkdir -p /gapped/AA
	mkdir -p /gapped/NT
	
### 5.1 Format paralog fasta files
#### **This example is from the output of 'paralog_investigator' of [HybPiper](https://github.com/mossmatters/HybPiper/wiki/Paralogs)**

To use the script from [Phylogenomic dataset construction respository](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/master/) an "@" symbol needs to be added in the fasta headers to identify paralogs of the same sample.

The fasta files *.paralogs.fasta after running 'paralog_investigator' has a format from [SPAdes](http://cab.spbu.ru/software/spades/) like

\>Rosa_woodsii  
AGTC...  
\>Lachemilla_pinnata.0 NODE_2_length_1787_cov_230.693567,Fragaria-gene15996_1557_01,4,519,78.47,(+),1,1537  
ACGT.....  
\>Alchemilla_colura.main NODE_3_length_1706_cov_62.896426,Fragaria-gene15996_1557_01,0,517,81.43,(-),1552,1  
ACCC....  
\>Alchemilla_colura.1 NODE_1_length_2101_cov_47.514174,Fragaria-gene15996_1557_01,0,519,79.11,(+),136,1687  
ACCG....  

##### To format the sequences run the loop below in the directory where the fasta files are locaded. Note: This will overwrite the fasta files.

an "@" symbol needs to be added in the fasta headers to identify paralogs of the same sample. To format the sequences run the loop below in the directory where the fasta files are locaded. Note: This will overwrite the fasta files from paralog investigator
run rename_paralogs.sh script in paralogs directory:

	for i in $(ls *.fasta); do 
	sed -i -E 's/(>.+)\.(.+)\s(.+)/\1@paralog_\2/' $i
	sed -i '/>/ {/@/! s/$/@unique/}' $i
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
In these steps, by distributing the genome outgroups in the hybpiper directory structure, those samples are incorporated when genes are gathered for aligment and named @main 
If other additional sequences are added to the fasta files at this stage (e.g. from reference genomes) make sure that those also have the "@" format.

### 5.2 align sequences using OMM_MACSE and build homolog trees
MACSE takes codon structure into account and runs efficiently with many samples
note: one issue that I ran into with the OMM_MACSE pipeline install, at least on the workstation I used, is that it has trouble with accessing directories other than the directory where it is installed. 
therefore, for large data sets I recommend installing the OMM_MACSE pipeline in the project directory, rather than your apps directory (otherwise you will have to copy all of your assemblies to your apps directory and mv the output alignment directories back to your project directory
run bash_MACSE.sh script from paralog_alignments folder to make individual bash files for each fasta gene file, this is the .sh script:

	for i in *paralogs.fasta; do 
	/{project_directory}/MACSE_V2_PIPELINES/OMM_MACSE/omm_macse_v10.02.sif 
	--in_seq_file $i 
	--out_dir ${i%%_FNA.fasta} 
	--out_file_prefix ${i%%_FNA.fasta}; 
	done > omm_macse_log.txt

Visually inspect alignments for artifacts! I removed artifacts manually in jalview. 
remove any genes that MACSE cant align (6128 is a problem for caryophyllaceae):
remove from the paralogs_unique.txt list and move the fasta file from the paralogs directory into another location (make a directory called "cant_align"?).

### 3.3 Replace "!" codon with gaps
If there are frame shifts in the alignments, MACSE will replace the shifted codons with "!" and will cause problems with RAxML or IQtree. 

	python /target_enrichment_orthology/scripts/remove_shifted_codons_from_macse.py /{project_directory}/hybpiper_output/orthology/ .AA.aln ./gapped/AA aa
	python /target_enrichment_orthology/scripts/remove_shifted_codons_from_macse.py /{project_directory}/hybpiper_output/orthology/ .NT.aln ./gapped/NT nt

### 3.4 trim alignments with Phyx
The output will files with extension .aln-cln

	python pxclsq_wrapper.py <alignment directory > <mininum column occupancy> <dna or aa>
	
	python /home/bencooper/Projects/target_enrichment_orthology/scripts/pxclsq_wrapper.py /media/data/cooper/PAFTOL/caryophyllaceae/hybpiper/paralogs/gapped/AA/ .90 aa
	python /home/bencooper/Projects/target_enrichment_orthology/scripts/pxclsq_wrapper.py /media/data/cooper/PAFTOL/caryophyllaceae/hybpiper/paralogs/gapped/NT/ .90 dna

phyx does wierd things to the fasta headers - 

### 3.5 Build homolog trees
infer trees with RAxML. This will infer trees with GTRGAMMA and 100 bootstrap replicates (different model and # of bs replicates can be modify in the script).

	python raxml_bs_wrapper.py <aln-cln files directory> <# of threads> <dna or aa>
	
Inspect your trees after making them 

### 3.6 Mask both mono- and paraphyletic tips that belong to the same taxon. 
Keep the tip that has the most un-ambiguous characters in the trimmed alignment. If you choose "n" for the parameter "mask_paraphyletic", it will only mask monophyletic tips. 

	python mask_tips_by_taxonID_transcripts.py <tre files directory> <aln-cln files directory> mask_paraphyletic <y or n>
	
Alternatively, to remove only monophyletic tips and keep the sequence with the shortest terminal branch length.

	python mask_tips_by_taxonID_genomes.py <tre files directory>
	

### 3.7 Trim spurious tips with [TreeShrink](https://github.com/uym2/TreeShrink)
if necessesary, otherwise procede to step 3.8

	python tree_shrink_wrapper.py <input directory> <tree file extension> <quantile>

It outputs the tips that were trimmed in the file .txt and the trimmed trees in the files .tt. You would need to test different quantiles to see which one fit better your data. Make sure you open the tree file to check whether TreeShrink removes your outgroups. 


##### These are you final homologs trees. If you want to make a final tree inference, write fasta files from trimmed trees. 

	python write_homolog_fasta_from_multiple_aln.py <original "@" formated fasta files directory> <trimmed trees directory> <fasta file extension> <tree file extension> <output directory>


## Step 3: Paralogy pruning to infer orthologs. Use one of the following


#### For details of orthology inference refer to Yang and Smith (2014) [doi:10.1093/molbev/msu245](https://doi.org/10.1093/molbev/msu245)



##### 1to1: filter homologs that contains only one sequence per species. No cutting is carried out.
	
	python filter_1to1_orthologs.py <final homologs directory> <tree file extension> <minimal number of taxa> <output directory>


##### MI: prune by maximum inclusion. Set OUTPUT_1to1_ORTHOLOGS to False if wish only to ouput orthologs that is not 1-to-1, for example, when 1-to-1 orthologs have already been analyzed in previous steps.

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



