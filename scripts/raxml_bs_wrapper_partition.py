"""
Input: a dir of cleaned alignments in fasta format and end with "-gb"
Output: trees estimated by raxml
"""

import os,sys
import subprocess
from seq import read_fasta_file

def raxml_bs(DIR,cleaned,num_cores,seqtype,replicates=200):
        assert cleaned.endswith(".aln-gb"),\
                "raxml infile "+cleaned+" not ends with .aln-gb"
        assert seqtype == "aa" or seqtype == "dna","Input data type: dna or aa"
        assert len(read_fasta_file(DIR+cleaned)) >= 4,\
                "less than 4 sequences in "+DIR+cleaned
        clusterID = cleaned.split(".")[0]
        tree = DIR+clusterID+".raxml_bs.part.tre"
        raw_tree = "RAxML_bipartitions."+cleaned
        model = "PROTCATAUTO" if seqtype == "aa" else "GTRGAMMA"
        partition = DIR+cleaned+".partition"
        if not os.path.exists(tree) and not os.path.exists(raw_tree):
                # raxml crashes if input file starts with .
                infasta = cleaned if DIR == "./" else DIR+cleaned
                
                cmd = ["raxml","-T",str(num_cores),\
                           "-f","a","-x","12345","-#",str(replicates),\
                           "-p","12345","-s",infasta,"-n",cleaned,"-m",model,"-q",partition]
                print " ".join(cmd)
                p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
                p.communicate()
                assert p.returncode == 0,"Error raxml"
        try:
                os.rename(raw_tree,tree)
                os.rename("RAxML_bootstrap."+cleaned, DIR+clusterID+".raxml_bs.part.trees")
                os.remove("RAxML_bestTree."+cleaned)
                os.remove("RAxML_info."+cleaned)
                os.remove("RAxML_log."+cleaned)
                os.remove("RAxML_parsimonyTree."+cleaned)
                os.remove("RAxML_result."+cleaned)
                os.remove("RAxML_bipartitionsBranchLabels."+cleaned)
                os.remove(DIR+cleaned+".reduced")
        except: pass # no need to worry about extra intermediate files
        return tree

def main(DIR,num_cores,seqtype):
        if DIR[-1] != "/": DIR += "/"
        filecount = 0
        for i in os.listdir(DIR):
                if i.endswith(".aln-gb"):
                        filecount += 1
                        raxml_bs(DIR,i,num_cores,seqtype)
        assert filecount > 0, "No file end with .aln-gb found in "+DIR


if __name__ == "__main__":
        if len(sys.argv) != 4:
                print "python raxml_bs_wrapper.py DIR number_cores dna/aa"
                print "make sure that the executable is named 'raxml' and is in the path"
                sys.exit(0)

        DIR,num_cores,seqtype  = sys.argv[1:]
        main(DIR,num_cores,seqtype)
