from Bio import AlignIO
import sys
import os
import getopt

# This small script takes a directory full of aligned fasta-files and trims those files: 
# 1. It removes all sequence parts in the beginning/end that aren't available for all sequences
# 2. It removes only keeps those parts of the aligments where sequences are avail. for each species
# 	 that are at least -l bases long (so it gets rid of small, wrongly aligned sequence-parts)
# 3. It trims the resulting alignments by -t bases. 
#
# The best: It keeps the reading frame! So you can use the resulting alignments for further, codon-sensitive, analysis
#
# Personal Use Case: I worked with a dataset that combines RNAseq-data with complete transcript-sequences of a couple of species. 
# As RNAseq-sequences tend to be much shorter and have quality problems, especially at the start/end of sequences, this led to 
# problems in down-stream (dN/dS) analysis. This script provided an easy way to get rid of misaligned sequence parts as well as
# of longer sequence-stretches that were only avail. for some species and due to RNAseq were missing for other species    

def read_alignment(infile):
	alignment = AlignIO.read(open(infile),"fasta")
	return alignment


def find_longest_start_gap(sequence):
	gap_counter = 1
	
	while gap_counter < len(sequence):
		if sequence.find("---"*gap_counter) == 0: 
			gap_counter += 1
		else:
			break
			
	gap_counter = gap_counter-1
	return gap_counter*3


def find_longest_end_gap(sequence):
	gap_counter = 1
	
	while gap_counter < len(sequence):
		if sequence[-3*gap_counter:] == "---"*gap_counter:
			gap_counter += 1
		else:
			break
	
	gap_counter = gap_counter-1
	
	return gap_counter*3


def row_iterator(alignment):
	start = 0
	end = 0
	
	for row in alignment:
		temp_start = find_longest_start_gap(str(row.seq))
		temp_end = find_longest_end_gap(str(row.seq))
		if temp_start > start:
			start = temp_start
		if temp_end > end:
			end = temp_end
	return [start,end]


def trim(path,infile):
	alignment = read_alignment(path+infile)
	delimitors = row_iterator(alignment)
	
	if delimitors[1] != 0:
		trimmed_alignment = alignment[:,delimitors[0]:-delimitors[1]]
	else:
		trimmed_alignment = alignment[:,delimitors[0]:]
	removed_bases(alignment,delimitors)	
	return trimmed_alignment


def removed_bases(alignment,delimitors):
	for row in alignment:
		print "Removed bases in front of "+row.id+": "+str(len(row.seq[:delimitors[0]])-row.seq[:delimitors[0]].count("-"))
		print "Removed bases in end of "+row.id+": "+str(len(row.seq[-delimitors[1]:])-row.seq[-delimitors[1]:].count("-"))


def remove_interim_gaps(alignment,cutoff):
	counter = 0
	hit_counter = 0
	new_alignment = []
	gap_counter = 0
	
	while counter < len(alignment[0]):
		base_pairing = alignment[:,counter]
		
		if base_pairing.find("-") == -1:
			hit_counter += 1
			gap_counter = 0
				
		elif gap_counter < 3:
			gap_counter += 1
			hit_counter += 1
		
		else:
			#CHECK FOR %3 in LINE 89 AND 97
			if hit_counter >= cutoff:
				if new_alignment == []:
					print "-- n=1 --"
					print "HIT:"
					print "Counter: "+str(counter)
					print "Hit Counter. "+str(hit_counter)
					print "Gap Counter: "+str(gap_counter)
					
					overhead = alignment[:,counter-hit_counter:counter-gap_counter].get_alignment_length()%3
					print alignment[:,counter-hit_counter+overhead:counter-gap_counter]
					print "Overhead: "+ str(overhead)
					new_alignment = alignment[:,counter-hit_counter+overhead:counter-gap_counter]
					
				else:
					print "-- n>1 --"
					print "HIT:"
					print "Counter: "+str(counter)
					print "Hit Counter. "+str(hit_counter)
					print "Gap Counter: "+str(gap_counter)
					
					overhead = alignment[:,counter-hit_counter:counter-gap_counter].get_alignment_length()%3
					print "Overhead: "+ str(overhead)
					new_alignment = new_alignment + alignment[:,counter-hit_counter+overhead:counter-gap_counter]
			hit_counter = 0
			gap_counter = 0
		counter += 1
		
	if hit_counter >= cutoff:
		if new_alignment == []:
			print "-- START & END --"
			print "HIT:"
			print "Counter: "+str(counter)
			print "Hit Counter. "+str(hit_counter)
			print "Gap Counter: "+str(gap_counter)
			
			overhead = alignment[:,counter-hit_counter:counter-gap_counter].get_alignment_length()%3
			print alignment[:,counter-hit_counter+overhead:counter-gap_counter]
			print "Overhead: "+ str(overhead)
			new_alignment = alignment[:,counter-hit_counter+overhead:counter-gap_counter]
		else:
			print "-- END --"
			print "HIT:"
			print "Counter: "+str(counter)
			print "Hit Counter. "+str(hit_counter)
			print "Gap Counter: "+str(gap_counter)
			overhead = alignment[:,counter-hit_counter:counter-gap_counter].get_alignment_length()%3
			print "Overhead: "+ str(overhead)
			new_alignment = new_alignment + alignment[:,counter-hit_counter+overhead:counter-gap_counter]
	hit_counter = 0
	
	return new_alignment

	
def alignment_writer(alignment_to_write,name,path,trim_length,file_type):
	if alignment_to_write != []:
		alignment_to_write = alignment_to_write[:,int(trim_length):-int(trim_length)]
		
		if alignment_to_write.get_alignment_length() >= 300:
			if alignment_to_write.get_alignment_length()%3 == 0:
				print "Length Alignment % 3: "+str(alignment_to_write.get_alignment_length()%3)
				AlignIO.write(alignment_to_write, path+name.replace(file_type,".msa.fasta"), "fasta")	
				AlignIO.write(alignment_to_write, path+name.replace(file_type,".msa.phylip"), "phylip")	
			else:
				print "ERROR "+name
		else:
			if alignment_to_write.get_alignment_length()%3 == 0:
				print "Length Alignment % 3: "+str(alignment_to_write.get_alignment_length()%3)
				AlignIO.write(alignment_to_write, path+name.replace(file_type,"_less_300bp.msa.fasta"), "fasta")	
				AlignIO.write(alignment_to_write, path+name.replace(file_type,"_less_300bp.msa.phylip"), "phylip")
			else:
				print "ERROR: "+name

def get_file_list(path,file_type):
        files = os.listdir(path)
        out_files = []
        for f in files:
                if f.find(file_type) != -1:
                        out_files.append(f)
        return out_files

def get_arguments(arguments):
	out_args = {}
	if len(arguments) == 1:
		print "use -h to get more help on the available options"
		print "e.g. ./trim_alignments.py -h"
	else:
		arguments.pop(0)
		optlist, arguments = getopt.getopt(arguments, 'hl:t:i:f:')
		for arg,opt in optlist:
			print arg
			print opt
			if arg == "-l":
				out_args["min_length_cutoff"] = opt
			elif arg == "-t":
				out_args["trim_length"] = opt
			elif arg == "-f":
				if opt[0] == ".":
					out_args["file_extension"] = opt
				else:
					out_args["file_extension"] = "."+opt
			elif arg == "-i":
				if opt[-1] != "/":
					out_args["infiles"] = opt+"/"
				else:
					out_args["infiles"] = opt
			elif arg == "-h":
				print "Welcome to the codon-alignment-trimmer!"
				print "You've got the following options: \n"
				print "-i, the folder with sequences that shall be processed (required)"
				print "-f, the file extension that shall be worked on (standard = fasta)"
				print "-l, minimum length that nucleotide-islands between gaps need to span (standard = 51)"
				print "-t, how many nucleotides shall be trimmed in front/at the end of the alignments (standard = 15)"
				print "-h, show this help message"
				break
		if out_args.has_key("trim_length") == False:
			out_args["trim_length"] = 15
		if out_args.has_key("min_length_cutoff") == False:
			out_args["min_length_cutoff"] = 51
		if out_args.has_key("file_extension") == False:
			out_args["file_extension"] = ".fasta"
	print out_args
	return out_args

def main():
	arguments = get_arguments(sys.argv)
	
	if arguments.has_key("infiles") == False:
		get_arguments(["-h"])
		exit()

	alignments = get_file_list(arguments["infiles"],arguments["file_extension"])
	for alignment in alignments:
		trimmed_alignment = trim(arguments["infiles"],alignment)
		new_alignment = remove_interim_gaps(trimmed_alignment,arguments["min_length_cutoff"])
		alignment_writer(new_alignment,alignment,sys.argv[1],arguments["trim_length"],arguments["file_extension"])
		
main()