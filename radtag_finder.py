#!/usr/bin/python

import sys
from Bio import SeqIO

# Read a tab-seperated file of restriction-enzyme-names and their target sequence
# e.g.
# EcoRI	GAATTC
# Returns hash. Key = enzyme-name, Value = enzyme-sequence

def tag_reader(infile):
	tag = {}
	handle = open(infile,"r")
	lines = handle.readlines()
	for line in lines:
		line_array = line.split("\t")
		tag[line_array[0]] = line_array[1].rstrip()
	return tag


# Read fasta-file, returns hash:
# key = sequence-identifier, value = sequence-information

def sequence_reader(infile):
	fasta = open(infile,"rU")
	sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
	return sequences


# iterates over all tags that were submitted, finds all restriction-sites
# in each sequence, prints them as tab-seperated list
# e.g. 
# Contig1	EcoR1	147

def tag_finder(tags,sequences):
	out_hash = {}
	for single_tag in tags:
		for single_sequence in sequences:
			finder = 0
			while finder != -1:
				finder = str(sequences[single_sequence].seq).find(tags[single_tag],finder)
				if finder != -1:
					print single_sequence+"\t"+single_tag+"\t"+str(finder)
					finder += 1

def main():
	try:
		tags = tag_reader(sys.argv[1])
		sequences = sequence_reader(sys.argv[2])
		tag_finder(tags,sequences)
	except: 
		print "######################################################################"
		print "You need to pass 2 parameters to use this script:"
		print "1. A tab-seperated list of the restriction-enzymes & their target sequence"
		print "2. A fasta-file of the sequences you want scanned"
		print ""
		print "Example: ./radtag_finder.py enzymes.tsv sequences.fasta > restriction_sites_sequences.tsv"

main()
