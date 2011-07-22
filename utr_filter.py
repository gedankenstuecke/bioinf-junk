#!/usr/bin/python

# This script reads a nucleotide sequence-fasta-file and the 
# output of ORFFinder to divide the nucleotide-sequences into
# the cds-regions and the pre- and post-cds UTRs. 

# The output can be found in utrs.fasta and cds.fasta
# USAGE: ./utr-filter.py some.fasta orffinder_results.tsv

from Bio import SeqIO
from Bio import Seq
import sys

# Read the FASTA-file, returns hash. Key = seq-identifier, value = Sequence-details

def sequence_reader(infile):
        fasta = open(infile,"rU")
        sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        return sequences


# Read the ORFPredictor-Output, returns hash. Key = seq-identifier, value = array.
# in the array:
# 1. if needed the reverse complement of the cds, to be able to search in sequences
# 2. the cds itself
# [maybe_reverse_complement, original_sequence]

def orf_reader(infile):
	orfs = {}
	handle = open(infile,"r")
	lines = handle.readlines()
	for line in lines:
		if line[0] != "#":
			line_array = line.split("\t")
			if int(line_array[1]) < 0:
				orfs[line_array[0]] = [Seq.reverse_complement(line_array[4]),line_array[4]]
			else:
				orfs[line_array[0]] = [line_array[4],line_array[4]]
	return orfs


# Search for UTRs in the original sequences

def utr_filter(sequences,orfs):
	utrs = {}
	cds = {}
	for orf in orfs:
		print orf
		sequence = str(sequences[orf].seq)
		if sequence[:sequence.find(orfs[orf][0])] != "":
			utrs[orf+"-UTR1"] = sequence[:sequence.find(orfs[orf][0])]		# find UTR that is in front of CDS
		try:
			if sequence[sequence.find(orfs[orf][0])+len(orfs[orf][0])] != "":
				utrs[orf+"-UTR2"] = sequence[sequence.find(orfs[orf][0])+len(orfs[orf][0]):]	# find UTR that follows the CDS
		except:
			pass

		cds[orf] = orfs[orf][1]		
	
	utr_handle = open("utrs.fasta","w")
	cds_handle = open("orfs.fasta","w")
	for single in utrs:
		print utrs[single]
		utr_handle.write(">"+single+"\n")
		utr_handle.write(utrs[single]+"\n")
	for single in cds:
		cds_handle.write(">"+single+"\n")
		cds_handle.write(cds[single]+"\n")

def main():
	try:
		sequences = sequence_reader(sys.argv[1])
		orfs = orf_reader(sys.argv[2])
		utr_filter(sequences,orfs)
	except:
		print "#################################################################################"
		print "You need to pass 2 parameters to use this script:"
		print "1. A FASTA-file of nucleotide-sequences that needs to be divided into UTRs & CDSs"
		print "2. The tab-seperated output of ORFFINDER, that is used to divide the sequences"
		print ""
		print "Example: ./utr_filter.py nucleotide_sequences.fasta orffilter_results.tsv"
