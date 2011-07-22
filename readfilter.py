#!/usr/bin/python

#### THIS SCRIPT FILTERS A HUGE FASTA-FILE ACCORDING TO A FILTER-LIST SUPPLIED ####
####								 									       ####
#### JUST PROVIDE 3 ARGUMENTS: a) THE FILTERING-LIST b) THE FASTA-FILE	       ####	
#### c) DO YOU WANT TO GET ALL READS OF THE FILTERING LIST? 				   ####
#### OR ALL READS _NOT_ ON THE LIST?										   ####	

from Bio import SeqIO
import sys



read_ids = open(sys.argv[1],"r")
readlist = read_ids.readlines()
reads = "".join(readlist)

fasta = open(sys.argv[2],"rU")
records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

out_file = open(sys.argv[2]+".filtered","w")
out_hash = {}

if sys.argv[3] == "include":
	for i in records:
		if reads.find(i) != -1:
			out_hash[i] = str(records[i].seq)
elif sys.argv[3] == "exclude":
	for i in records:
		if reads.find(i) == -1:
			out_hash[i] = str(records[i].seq)
else:
	print "specify job as 'include' or 'exclude'"
	print "'include' keeps only sequences specified in the list-file"
	print " while 'exclude' keeps only sequences NOT specified in the list-file"

for i in out_hash:
	out_file.write(">"+i+"\n")
	out_file.write(out_hash[i])
	out_file.write("\n")

