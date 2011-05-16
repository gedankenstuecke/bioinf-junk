#!/usr/bin/python

#### THIS SCRIPT FILTERS A HUGE FASTA-FILE ACCORDING TO A FILTER-LIST SUPPLIED ####
#### 																	       ####
#### JUST PROVIDE 3 ARGUMENTS: a) THE FILTERING-LIST b) THE FASTA-FILE	       ####
#### i.e. ./readfilter.py input.fasta filtering.txt include					   ####	

from Bio import SeqIO
import sys



read_ids = open(sys.argv[1],"r")
readlist = read_ids.readlines()
reads = {}

for i in readlist:
	x = i.split("\t")
	est = x[0]
	est = est.replace("\"","")
	contig = x[1]
	contig = contig.replace("\"","")
	contig = contig.rstrip("\n")
	if contig in reads:
		reads[contig].append(est)
	else:
		reads[contig] = [est]
print reads

contigfile = open(sys.argv[2],"rU")
contigsequence = SeqIO.to_dict(SeqIO.parse(contigfile, "fasta"))

estfile = open(sys.argv[3],"rU")
estsequence = SeqIO.to_dict(SeqIO.parse(estfile, "fasta"))

for i in reads:
	outcontig = open(i+".fasta","w")
	outreads = open(i+".fasta.reads","w")
	outcontig.write(">"+i+"\n")
	outcontig.write(str(contigsequence[i].seq)+"\n")
	for j in reads[i]:
		outreads.write(">"+j+"\n")
		outreads.write(str(estsequence[j].seq)+"\n")
	outcontig.close()
	outreads.close()

#fasta = open(sys.argv[2],"rU")
#records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

#out_file = open(sys.argv[2]+".filtered","w")
#out_hash = {}

#for i in readlist:
#	i = i.rstrip()
#	out_hash[i] = str(records[i].seq)

#for i in records:
#	if reads.find(i) == -1:
#		out_hash[i] = str(records[i].seq)

#for i in out_hash:
#	out_file.write(">"+i+"\n")
#	out_file.write(out_hash[i])
#	out_file.write("\n")

