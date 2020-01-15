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
