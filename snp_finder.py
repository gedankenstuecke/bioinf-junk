#!/usr/bin/python
# This script allows to call SNPs out of a bam-alignment-file by looking at 
# the depth of reads at a given position in a chromosome/contig and the 
# allele-frequency of the bases at this position

import pysam
import sys

# open the BAM-file 

def open_sam(infile):
	samfile = pysam.Samfile(infile, "rb")
	return samfile	


# get the SNPs according to users parameters

def get_snps(samfile,min_coverage,min_frequency):
	chromosome_name = ""
	out_snps = []
	
	# iterate over all contigs & all positions	
	
	for pileupcolumn in samfile.pileup(): 	
		
		if chromosome_name != samfile.getrname(pileupcolumn.tid):
			print "Scanning %s" % (samfile.getrname(pileupcolumn.tid))
			chromosome_name = samfile.getrname(pileupcolumn.tid)
		else:
			pass
		
		position = pileupcolumn.pos			# get current position
		reads_at_position = pileupcolumn.n	# get amount of reads at this position
		
		if reads_at_position >= int(min_coverage):	#check that minimum coverage is met
			potential_snp = snp_counter(pileupcolumn,position,min_frequency,chromosome_name,reads_at_position)
			if len(potential_snp.genotypes) > 1:		# only keep SNPs with at least 2 genotypes 
				out_snps.append(potential_snp)
	return out_snps


# iterates over all reads in a column to calculate nucleotide-frequencies 
# returns one snp-class-item

def snp_counter(pileupcolumn,position,min_frequency,chromosome_name,reads_at_position):
	bases = {"A":0,"T":0,"C":0,"G":0}
	genotypes = {}
	
	# iterate over all reads 
	# add to base-cound if a standard-base
	
	for pileupread in pileupcolumn.pileups:		
		if "ATGC".find(pileupread.alignment.seq[pileupread.qpos]) != -1:
			bases[pileupread.alignment.seq[pileupread.qpos]] += 1	
	
	#check that the minimum allele frequency is met
	# if so: add to snp-class
	
	for base in bases:
		if float(bases[base])/reads_at_position > float(min_frequency):
			genotypes[base] = float(bases[base])/reads_at_position
	snp_return = snp(chromosome_name,position+1,genotypes)
	return snp_return			


# print SNPs as a tsv-file:

def snp_writer(snp_array):
	out_handle = open("snps.tsv","w")
	
	for single_snp in snp_array:
		snp_line = single_snp.chromosome+"\t"+str(single_snp.position)
		for genotype in single_snp.genotypes:
			snp_line = snp_line + "\t" + genotype + "\t" + str(single_snp.genotypes[genotype])
		snp_line = snp_line + "\n"
		out_handle.write(snp_line)


def main():
	try:
		samfile = open_sam(sys.argv[1])
		snps = get_snps(samfile,sys.argv[2],sys.argv[3])
		snp_writer(snps)
	except:
		print "######################################################################"
		print "You need to pass 3 parameters to use this script:"
		print "1. A sorted and indexed BAM-file (produced by bwa/bowtie and samtools)"
		print "2. The minimum depth of reads that needs to be achieved to call a SNP"
		print "3. The minimum allele frequency of a base to be counted as a SNP"
		print ""
		print "Example: ./snp_finder.py input.bam 10 0.2"

	
class snp():
	def __init__(self,chromosome,position,genotypes):
		self.chromosome = chromosome
		self.position = position
		self.genotypes = genotypes

main()
