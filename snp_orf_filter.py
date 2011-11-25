import sys
import getopt
from Bio import SeqIO
from Bio.Seq import Seq

## This script takes 3 input-files to get a list of SNPs that occur inside ORFs, which where predicted using ORFFinder
## Those SNPs get checked to see if they are synonymous/non-synonymous. Results will be split in two files for the categories.
## 
## Files the user needs to deliver: 
## 1. Original Fasta-file of the contigs/chromosomes
## 2. ORFFinder-result-tsv-file
## 3. SNP-report-file 

# Get all Parameters
def get_arguments(arguments):
	"""Parse all the parameters we've got"""
	out_args = {}
	if len(arguments) == 1:
		print "Run with -h to get help"
		exit()
	else:
		arguments.pop(0)
		optlist, arguments = getopt.getopt(arguments, 'hi:s:o:f:')
		for arg,opt in optlist:
			if arg == "-h":
				print "Parameters:"
				print "-i: Tabular ORF-File (required)"
				print "-s: SNP-Report-File (required)"
				print "-o: Output-File (standard: snps_filtered.tsv)"
				print "-f: untrimmed Fasta-file with UTR & ORF (required)"
				exit()
			elif arg == "-i":
				out_args["orfs"] = opt
			elif arg == "-s":
				out_args["snps"] = opt
			elif arg == "-o":
				out_args["outfile"] = opt
			elif arg == "-f":
				out_args["original_contigs"] = opt
		if out_args.has_key("orfs") == False:
			print "Please provide an ORF-file via -i"
			exit()
		if out_args.has_key("snps") == False:
			print "Please provide a SNP-file via -s"
			exit()
		if out_args.has_key("original_contigs") == False:
			print "Please provide a sequence-fasta-file via -f"
			exit()
		return out_args

# SNP-Class that saves contig, position & full line of original file
class snp():
	def __init__(self,chromosome,position,full_line,minor_allele):
		self.chromosome = chromosome
		self.position = position
		self.full_line = full_line
		self.minor_allele = minor_allele
	


# READ SNP-File	
def read_snps(snp_file):
	snp_array = []
	handle = open(snp_file,"r")
	snps = handle.readlines()
	snps.pop(0)
	handle.close()
	for line in snps: 
		temp_array = line.rstrip().split("\t")
		snp_array.append(snp(temp_array[1],int(temp_array[2]),temp_array,temp_array[6]))
	return snp_array		

# ORF-Class that saves name, start & stop
class orf():
	def __init__(self,name,start,stop,orientation):
		self.name = name
		self.start = start
		self.stop = stop
		self.orientation = orientation
	


# READ ORF-File
def read_orfs(orf_file):
	orf_hash = {}
	handle = open(orf_file,"r")
	orf_lines = handle.readlines()
	orf_lines.pop(0)
	handle.close()
	for line in orf_lines:
		temp_array = line.split("\t")
		if int(temp_array[2]) > int(temp_array[3]):
			orf_hash[temp_array[0]] = orf(temp_array[0], int(temp_array[3]), int(temp_array[2]), "-")
		else:
			orf_hash[temp_array[0]] = orf(temp_array[0], int(temp_array[2]), int(temp_array[3]), "+")
	return orf_hash

# Find all SNPs that are located inside the ORFs
def snp_in_orf(unfiltered_snps, orfs, sequences):
	syn_snps = []
	non_syn_snps = []
	for single_snp in unfiltered_snps:
		if orfs.has_key(single_snp.chromosome):
			if single_snp.position > orfs[single_snp.chromosome].start:
				if single_snp.position < orfs[single_snp.chromosome].stop: 
					print "SNP inside of ORF, start checking for synonymous/non-synonymous polymorphism"
					syn_tuple = syn_check(single_snp,orfs[single_snp.chromosome],sequences)
					print "Checked, SNP is %s" %(syn_tuple[0])
					if syn_tuple[0] == "syn":
						syn_snps.append(syn_tuple[1])
					else:
						non_syn_snps.append(syn_tuple[1])
	print "Got all SNPs that are located inside ORFs"
	out_hash = {"non-syn":non_syn_snps,"syn":syn_snps}
	return out_hash


# Check if SNPs are syn/non-syn and return tuple (syn/nonsyn, snp). snp has 2 more fields in full_line: original codon, polymorph codon 
def syn_check(single_snp,single_orf,sequences):
	if single_orf.orientation == "+":
		print "HIT: "+single_snp.chromosome
		codon_position = (single_snp.position - single_orf.start)%3 # (total position of snp - stuff that gets cut away) = relative position of snp in codon-terms 
		if codon_position == 1:
			print "FIRST POSITION"
			original_codon = str(sequences[single_snp.chromosome].seq[single_snp.position:single_snp.position+3])
			minor_codon = str(single_snp.minor_allele+original_codon[1:])
			print "ORIGINAL CODON: "+str(original_codon)
			print "MINOR CODON: "+str(minor_codon)
		elif codon_position == 2: 
			print "SECOND POSITION"
			original_codon = str(sequences[single_snp.chromosome].seq[single_snp.position-1:single_snp.position+2])
			minor_codon = str(sequences[single_snp.chromosome].seq[single_snp.position-1]+single_snp.minor_allele+sequences[single_snp.chromosome].seq[single_snp.position+1])
			print "ORIGINAL CODON: "+str(original_codon)
			print "MINOR CODON: "+ minor_codon
		else:
			print "THIRD POSITION"
			original_codon = str(sequences[single_snp.chromosome].seq[single_snp.position-2:single_snp.position+1])
			minor_codon = str(original_codon[:2]+single_snp.minor_allele)
			print "ORIGINAL CODON: "+str(original_codon)
			print "MINOR CODON: "+minor_codon
	else:
		# think about how to get position for backwards-running orfs, somehow think that it should be: single_orf.stop (this is base #1) - 
		codon_position = (single_orf.stop - single_snp.position)%3
		if codon_position == 1: 
			print "REVERSE #1"
			original_codon = str(sequences[single_snp.chromosome].seq[single_snp.position]+sequences[single_snp.chromosome].seq[single_snp.position-1]+sequences[single_snp.chromosome].seq[single_snp.position-2])
			minor_codon = single_snp.minor_allele+sequences[single_snp.chromosome].seq[single_snp.position-1]+sequences[single_snp.chromosome].seq[single_snp.position-2]
			print original_codon
			print minor_codon
		elif codon_position == 2:
			print "REVERSE #2" 
			original_codon = str(sequences[single_snp.chromosome].seq[single_snp.position+1]+sequences[single_snp.chromosome].seq[single_snp.position]+sequences[single_snp.chromosome].seq[single_snp.position-1])
			minor_codon = str(sequences[single_snp.chromosome].seq[single_snp.position+1]+single_snp.minor_allele+sequences[single_snp.chromosome].seq[single_snp.position-1])
			print original_codon
			print minor_codon
		else:
			print "REVERSE #3"
			original_codon = str(sequences[single_snp.chromosome].seq[single_snp.position+2]+sequences[single_snp.chromosome].seq[single_snp.position+1]+sequences[single_snp.chromosome].seq[single_snp.position])
			minor_codon = str(sequences[single_snp.chromosome].seq[single_snp.position+2]+sequences[single_snp.chromosome].seq[single_snp.position+1]+single_snp.minor_allele)
			print original_codon
			print minor_codon
	
	original_codon_translated = Seq(original_codon).translate()
	minor_codon_translated = Seq(minor_codon).translate()
	print original_codon_translated
	print minor_codon_translated
	single_snp.full_line.append(original_codon)
	single_snp.full_line.append(minor_codon)
	if str(original_codon_translated) == str(minor_codon_translated):
		return ("syn",single_snp)
	else:
		return ("non-syn",single_snp)


# READ FASTA FILE
def sequence_reader(infile):
	fasta = open(infile,"rU")
	sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
	return sequences

# WRITE SNPs
def write_snps(snp_hash,outfile_name):
	print "Writing SNP-files"
	handle = open(outfile_name+"_synonymous.txt","w")
	syn_snps = snp_hash["syn"]
	for single_snp in syn_snps:
		handle.write("\t".join(single_snp.full_line)+"\n")
	handle.close()
	handle = open(outfile_name+"_non_synonymous.txt","w")
	non_syn_snps = snp_hash["non-syn"]
	for single_snp in non_syn_snps:
		handle.write("\t".join(single_snp.full_line)+"\n")
	handle.close()
	print "Written files"

# Make all that stuff run in a "sane" fashion
def main():
	"""Run all the stuff"""
	arguments = get_arguments(sys.argv)
	unfiltered_snps = read_snps(arguments["snps"]) 
	sequences = sequence_reader(arguments["original_contigs"])
	print unfiltered_snps[0].chromosome
	print unfiltered_snps[0].position
	print unfiltered_snps[0].full_line
	orfs = read_orfs(arguments["orfs"])
	snp_hash = snp_in_orf(unfiltered_snps,orfs,sequences)
	write_snps(snp_hash, arguments["outfile"])

main()