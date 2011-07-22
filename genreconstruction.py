#!/usr/bin/python

from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
import sys
import re

# TRY TO OPEN THE INPUT FILE #

def input():
	try:
		infile = open(sys.argv[1])
		fastafile = open(sys.argv[2],"rU")
		fasta = SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))
		return [infile,fasta]
		
	except: 
		print "This script takes a standard NCBI-Blast XML-file and the corresponding FASTA-file."
		print "It filters reads according to the evalue cutoff (1.0e-15) and renames the sequences."
		print "The renaming is done by taking the input-name and the blast-result-name."
		print "The sequences are shortened to the length of the BLAST-hit"
		print ""
		print "Just start the script giving a single XML-File and a FASTA-file"
		print "i.e. \"./volatiltyfilter.py blastout.xml input.fasta\""


def parser(handler,fasta):
	blast_records = NCBIXML.parse(handler)
	outhash = {}
	print ""
	print"FILTERING READS..."
	print""
		
	for i in blast_records:
		hitlist = []
		hitlist2 = []
		read = i.query	
		for alignment in i.alignments:
			aligntitle = alignment.title
			for hsp in alignment.hsps:
				expect = hsp.expect
				identity = hsp.identities
				start = hsp.query_start
				query = len(hsp.query*3)
				seq = str(fasta[read].seq)
				seq = seq[start-1:start+query-1]
				
				if expect < 1.0e-15:
					if hitlist != []:
						for i in range(0,len(hitlist)):
														
							# IF START OF CURRENT QUERY IS BETWEEN START AND STOP OF CURRENT HITLIST  
							# AND THE CURRENT QUERY IS LONGER AS THE CURRENT HIT REPLACE THE HITLIST WITH QUERY
							
							if hitlist[i][4] >= start:
								if start >= hitlist[i][3]:
									if query > hitlist[i][5]:
										hitlist[i] = [read,aligntitle,seq,start,start+query-1,query]
							#			print "REPLACE 1:" +  read + " " + str(start) + " " + str(start+query-1)
										
							# IF STOP OF CURRENT QUERY IS BETWEETN START AND STOP OF CURRENT HITLIST
							# AND THE CURRENT QUERY IS LONGER AS THE CURRENT HIT REPLACE THE HITLIST WITH QUERY			
										
							elif hitlist[i][4] >= start+query-1:
								if start+query-1 >= hitlist[i][3]:
									if query > hitlist[i][5]:
										hitlist[i] = [read,aligntitle,seq,start,start+query-1,query]
							#			print "REPLACE 2:" + read + " " + str(start) + " " + str(start+query-1)
										
							# FIND OUT IF NEW QUERY IS NOT ALREADY IN HITLIST AND REALLY IS A NEW HIT
							# CHECK ALL HITLISTS AND IF SUM OF NEGATIVE HITS == LEN(HITLIST) -> ADD IT
							
							else:
								hitcount = 0
								for j in range(0,len(hitlist)):
									if start+query-1 < hitlist[j][3]:
										if start < hitlist[j][3]:
											hitcount += 1
									elif start+query-1 > hitlist[j][4]:
										if start > hitlist[j][4]:
											hitcount += 1
								if hitcount == len(hitlist):
									hitlist.append([read,aligntitle,seq,start,start+query-1,query])
							#		print "APPEND HIT: "+ read + " " + str(start) + " " + str(start+query-1)
									
					else:
						hitlist.append([read,aligntitle,seq,start,start+query-1,query])
						#print "NEW: " + read + " " + str(start) + " " + str(start+query-1)
						
		if hitlist != []:
			hitlist = hitlist + hitlist2
			for i in hitlist:
				outhash[i[0]+" "+i[1]] = i[2]
	
		print "DONE FILTERING"
		print ""
			
	return outhash
	

def writer(input):
	outfile = open(sys.argv[2]+"filtered_queries.fasta","w")
	ginames = []
	print "WRITING FILTERED READS"
	for i in input:
		outfile.write(">"+i+"\n")
		outfile.write(input[i])
		outfile.write("\n")
		gistart = i.find("gi|")
		gistop = i.find("|",gistart+3)
		ginames.append(i[gistart+3:gistop])
	print ""
	print "DONE WRITING READS"
	return ginames

def gifinder(ginames):
	print ""
	print "START FILTERING HITS"
	outthing = {}
	Entrez.email = "bgreshake@googlemail.com"
	oldname = ""
	outfile2 = open(sys.argv[2]+"filtered_hits.fasta","w")
	for i in ginames:
		handle = Entrez.efetch(db="protein", id = i, rettype = "gb")
		record = SeqIO.read(handle,"genbank")
		for j in record.features:
			if j.type == "CDS":
				identifier = j.qualifiers["coded_by"][0]
				identexpression = re.compile("[a-zA-Z0-9_.]+\:")
				startexpression = re.compile("[<>\d]+\.\.")
				stopexpression = re.compile("\.\.[<>\d]\d*")				
				sequence = ""				
				if identifier.find("complement") != -1:
					if identifier.find("join") != -1: 
						identifier = identifier.split(",")
						for k in identifier: 
							ident = identexpression.search(k)
							name = ident.group(0)	
							name = name[:-1]
							restart = startexpression.search(k)
							start = restart.group(0)
							start = re.sub("\D","",start)
							restop = stopexpression.search(k)
							stop = restop.group(0)
							stop = re.sub("\D","",stop)
							if oldname != name:
								new_handle = Entrez.efetch(db="nucleotide", id = name, rettype = "fasta")
								new_record = SeqIO.read(new_handle,"fasta")
							tempsequence = str(new_record.seq)
							sequence = sequence + tempsequence[int(start)-1:int(stop)-1]
							oldname = name
						bioseq = Seq(sequence)
						sequence = str(bioseq.complement())				
						outfile2.write(">" + name)
						outfile2.write(sequence)
					else: 
						ident = identexpression.search(identifier)
						name = ident.group(0)
						name = name[:-1]
						restart = startexpression.search(identifier)
						start = restart.group(0)
						start = re.sub("\D","",start)
						restop = stopexpression.search(identifier)
						stop = restop.group(0)
						stop = re.sub("\D","",stop)
						if oldname != name:
							new_handle = Entrez.efetch(db="nucleotide", id = name, rettype = "fasta")
							new_record = SeqIO.read(new_handle,"fasta")
						tempsequence = str(new_record.seq)
						sequence = sequence + tempsequence[int(start)-1:int(stop)-1]
						oldname = name
						bioseq = Seq(sequence)
						sequence = str(bioseq.complement())
						outfile2.write(">" + name)		
						outfile2.write(sequence)
				else:
					if identifier.find("join") != -1:
                                                identifier = identifier.split(",")
                                                for k in identifier:
                                                        ident = identexpression.search(k)
                                                        name = ident.group(0)
                                                        name = name[:-1]
                                                        restart = startexpression.search(k)
                                                        start = restart.group(0)
                                                        start = re.sub("\D","",start)
                                                        restop = stopexpression.search(k)
                                                        stop = restop.group(0)
                                                        stop = re.sub("\D","",stop)                                               
                                                        if oldname != name:
                                                                new_handle = Entrez.efetch(db="nucleotide", id = name, rettype = "fasta")
                                                                new_record = SeqIO.read(new_handle,"fasta")
                                                        tempsequence = str(new_record.seq)
                                                        sequence = sequence + tempsequence[int(start)-1:int(stop)-1]
                                                        oldname = name
						outfile2.write(">" + name)
                                                outfile2.write(sequence)
					else:	
                                                ident = identexpression.search(identifier)
                                                name = ident.group(0)
                                                name = name[:-1]
                                                restart = startexpression.search(identifier)
                                                start = restart.group(0)
                                                start = re.sub("\D","",start)
                                                restop = stopexpression.search(identifier)
                                                stop = restop.group(0)
                                                stop = re.sub("\D","",stop)
                        	                if oldname != name:
							 new_handle = Entrez.efetch(db="nucleotide", id = name, rettype = "fasta")
                                                         new_record = SeqIO.read(new_handle,"fasta")
                                                tempsequence = str(new_record.seq)
                                                sequence = sequence + tempsequence[int(start)-1:int(stop)-1]
                                                oldname = name
						outfile2.write(">" + name)
                                                outfile2.write(sequence)
	print "DONE FILTERING & WRITING HITS"
	print ""
	print "HAVE A NICE DAY"

handler = input()

if handler != None:
	outlist = parser(handler[0],handler[1])
	output = writer(outlist)
	gifinder(output)	
