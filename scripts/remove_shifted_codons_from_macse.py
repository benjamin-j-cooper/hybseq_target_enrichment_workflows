###This script replaces codons and aa with frame shift from macse aligments with gaps.


from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os,sys


def codon(inDIR,file_ending,outDIR):

	if os.path.isabs(inDIR) == False: inDIR = os.path.abspath(inDIR)
	if inDIR[-1] != "/": inDIR += "/"
	if os.path.isabs(outDIR) == False: outDIR = os.path.abspath(outDIR)
	if outDIR[-1] != "/": outDIR += "/"
	if outDIR == ".": outDIR = os.getcwd()
	
	
	filecount = 0

	codon_shift_array = ["!!A", "!!C", "!!G", "!!T", "A!!", "C!!", "G!!", "T!!", "!A!", "!C!", "!G!", "!T!", "AC!", "AG!", "AT!", "AA!", "CA!", "CG!", "CT!", "CC!", "GA!", "GC!", "GT!", "GG!", "TA!", "TC!", "TG!", "TT!", "!AC", "!AG", "!AT", "!AA", "!CA", "!CG", "!CT", "!CC", "!GA", "!GC", "!GT", "!GG", "!TA", "!TC", "!TG", "!TT", "A!C", "A!G", "A!T", "A!A", "C!A", "C!G", "C!T", "C!C", "G!A", "G!C", "G!T", "G!G", "T!A", "T!C", "T!G", "T!T",
					"!!a", "!!c", "!!g", "!!t", "a!!", "c!!", "g!!", "t!!", "!a!", "!c!", "!g!", "!t!", "ac!", "ag!", "at!", "aa!", "ca!", "cg!", "ct!", "cc!", "ga!", "gc!", "gt!", "gg!", "ta!", "tc!", "tg!", "tt!", "!ac", "!ag", "!at", "!aa", "!ca", "!cg", "!ct", "!cc", "!ga", "!gc", "!gt", "!gg", "!ta", "!tc", "!tg", "!tt", "a!c", "a!g", "a!t", "a!a", "c!a", "c!g", "c!t", "c!c", "g!a", "g!c", "g!t", "g!g", "t!a", "t!c", "t!g", "t!t"]

	for i in os.listdir(inDIR):
		if i.endswith(file_ending):
			print i
			filecount += 1
						
			basename = os.path.splitext(i)[0]
			new_aln = outDIR+basename+".fs.aln"
			file=open(new_aln,'w')
			aln = SeqIO.parse(inDIR+i, "fasta")
			for record in aln:

				if "!" in record.seq:
					tempRecordSeq = list(record.seq)
					for index in range(0, len(record.seq), 3):
						codon = record.seq[index:index+3]
						if codon in codon_shift_array:
							tempRecordSeq[index:index+3] = "---"
			        	record.seq = Seq("".join(tempRecordSeq))
				else:
					record.seq == record.seq
					
				id =record.id
				seq=str(record.seq)
				file.write(">" + id + "\n" + seq + "\n")

			file.close()
			
	assert filecount > 0, \
		"No file end with "+file_ending+" found in "+inDIR
			
			


def aa(inDIR, file_ending, outDIR):

	if os.path.isabs(inDIR) == False: inDIR = os.path.abspath(inDIR)
	if inDIR[-1] != "/": inDIR += "/"
	if os.path.isabs(outDIR) == False: outDIR = os.path.abspath(outDIR)
	if outDIR[-1] != "/": outDIR += "/"
	if outDIR == ".": outDIR = os.getcwd()
	
	
	filecount = 0

	aa_shiflt_array = ["!", "*"]
	
	for i in os.listdir(inDIR):
		if i.endswith(file_ending):
			print i
			filecount += 1	
								
			basename = os.path.splitext(i)[0]
			new_aln = outDIR+basename+".fs.aln"
			aln = SeqIO.parse(inDIR+i, "fasta")
			file=open(new_aln,'w')
			for record in aln:
				if "!" or "*" in record.seq:
					tempRecordSeq = list(record.seq)
					for index in range(0, len(record.seq), 1):
						aa = record.seq[index:index+1]
						if aa in aa_shiflt_array:
							tempRecordSeq[index:index+1] = "-"
						record.seq = Seq("".join(tempRecordSeq))
				else:
					record.seq == record.seq
					
				id =record.id
				seq=str(record.seq)
				file.write(">" + id + "\n" + seq + "\n")

			file.close()

			
	assert filecount > 0, \
		"No file end with "+file_ending+" found in "+inDIR
			

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "Usage:"
		print "python remove_shifted_codons_from_macse.py inDIR, file_ending, outDIR, alingment_type [nt or aa]"
	else:	
		assert sys.argv[4] == "nt" or sys.argv[4] == "aa"
		"alignment type has to be either nt or aa"
		if sys.argv[4] == "nt":
			codon(inDIR=sys.argv[1],file_ending=sys.argv[2],outDIR=sys.argv[3])

		else:
			aa(inDIR=sys.argv[1],file_ending=sys.argv[2],outDIR=sys.argv[3])
	sys.exit(0)		