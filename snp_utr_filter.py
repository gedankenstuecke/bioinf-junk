#!/usr/bin/python

# This script reads a nucleotide sequence-fasta-file, the output 
# of ORFFinder and a SNP-file to divide SNPs into
# SNPS in the coding region and the UTRs 

# The output can be found in snps_cds.tsv, snps_utr.tsv and snps_na.tsv
# USAGE: ./utr-filter.py some.fasta orffinder_results.tsv snps.tsv

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
	cds = {}
	for orf in orfs:
		sequence = str(sequences[orf].seq)
		start = sequence.find(orfs[orf][0]) + 1
		end = sequence.find(orfs[orf][0]) + 1 + len(orfs[orf][0])			
		cds[orf] = [start,end]
		
		print orf
		print start
		print end
	return cds
	
def snp_reader(infile):
		handle = open(infile,"r")
		lines = handle.readlines()
		return lines

def snp_divider(cds,lines):
	utr_snps = []
	cds_snps = []
	na_snps = []
	for line in lines:
		line_array = line.split("\t")
		if cds.has_key(line_array[0]):
			if int(line_array[1]) >= int(cds[line_array[0]][0]) and int(line_array[1]) <= int(cds[line_array[0]][1]):
				cds_snps.append(line_array)
			else:
				utr_snps.append(line_array)
		else:
			na_snps.append(line_array)
	
	utr_handle = open("snps_utr.tsv","w")
	for line_array in utr_snps:
		utr_handle.write("\t".join(line_array))
	cds_handle = open("snps_cds.tsv","w")
	for line_array in cds_snps:
		cds_handle.write("\t".join(line_array))
	na_handle = open("snps_na.tsv","w")
	for line_array in na_snps:
		na_handle.write("\t".join(line_array)) 

def main():
	try:
		sequences = sequence_reader(sys.argv[1])
		orfs = orf_reader(sys.argv[2])
		cds = utr_filter(sequences,orfs)
		snps = snp_reader(sys.argv[3])

	except:
		print "#################################################################################"
		print "You need to pass 2 parameters to use this script:"
		print "1. A FASTA-file of reference nucleotide-sequences"
		print "2. The tab-seperated output of ORFFINDER, that is used to divide the SNPs"
		print "3. The SNPs.tsv that shall be divided"
		print ""
		print "Example: ./utr_filter.py nucleotide_sequences.fasta orffilter_results.tsv snps.tsv"
	snp_divider(cds, snps)
main()