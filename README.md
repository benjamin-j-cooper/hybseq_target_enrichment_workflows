# Target enrichment orthology

#### This repository contains scripts for orthology in inference from target enrichment data (e.g. Hyb-Seq) and relies mainly in scripts from [Phylogenomic dataset construction respository](https://bitbucket.org/yangya/phylogenomic_dataset_construction/src/master/) plus some additional ones.

#### If using this scrips from this repository you must cite:

Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution. doi: 10.1093/molbev/msu245

Morales-Briones, D.F., G. Kadereit, D.T. Tefarikis, M.J. Moore, S.A. Smith, S.F. Brockington, A.Timoneda, W.C. Yim, J.C. Cushman, Y. Yang. 2019. Disentangling Sources of Gene Tree Discordance in Phylotranscriptomic Datasets: A Case Study from Amaranthaceae s.l. bioRxiv 794370.



## Dependencies needed to run the scripts. 

[TreeShrink](https://github.com/uym2/TreeShrink) It works now with Version 1.3.2 (older versions won't work)

[RAxML](https://github.com/stamatak/standard-RAxML) Version 8.2.11  (newer versions should work)

[Phyx](https://github.com/FePhyFoFum/phyx)

[mafft](https://mafft.cbrc.jp/alignment/software/) Version 7.307 (newer versions should work)

[macse](https://bioweb.supagro.inra.fr/macse/index.php?menu=releases) Version 2.0.3 (newer versions should work)



### Step 1: Format paralog fasta files

#### **This example if for the output of 'paralog_investigator' of [HybPiper](https://github.com/mossmatters/HybPiper/wiki/Paralogs)**

To use the script from [Phylogenomic dataset construction respository](https://bitbucket.org/yangya/phylogenomic_dataset_construction/src/master/) an "@" symbol needs to be added to identify paralogs of the same sample.

The fasta files *.paralogs.fasta after running 'paralog_investigator' have a format from [SPAdes](http://cab.spbu.ru/software/spades/) like:

\>Rosa_woodsii
AGTC...
\>Lachemilla_pinnata.0 NODE_2_length_1787_cov_230.693567,Fragaria-gene15996_1557_01,4,519,78.47,(+),1,1537
ACGT.....
\>Alchemilla_colura.main NODE_3_length_1706_cov_62.896426,Fragaria-gene15996_1557_01,0,517,81.43,(-),1552,1
ACCC....
\>Alchemilla_colura.1 NODE_1_length_2101_cov_47.514174,Fragaria-gene15996_1557_01,0,519,79.11,(+),136,1687
ACCGG....

To format the sequences run the next loop in the folder with the fasta files:

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
ACCGG....

For details of SPAdpes contigs see [HybPiper](https://github.com/mossmatters/HybPiper/wiki/Paralogs)


If other additional sequences are added to the fasta files (e.g. from reference genomes) make sure that those also have the "@" format.

I added reference sequences of <em>Fragaria vesca</em> as Fragaria_vesca@genome
