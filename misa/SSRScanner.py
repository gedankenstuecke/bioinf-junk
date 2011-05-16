#!/usr/bin/python
from Bio import SeqIO
import sys
import os

## START THE MAIN PROGRAM ##

def main():
	if len(sys.argv) != 2:
		print "Please try --help or -h for help"


	elif sys.argv[1] == "--help" or sys.argv[1] == "-h" or len(sys.argv) < 2:
		print ""
		print "###############################################"
		print "This tool finds microsatellites using MISA"
		print "The settings for misa can be found in misa.ini"
		print "The export is done using a csv-file"
		print ""
		print "Just pass over a single FASTA-file" 
		print "###############################################"
		print ""

	else:

		try:
			descriptorunifier(sys.argv[1]) 
			os.system("./misa.pl "+sys.argv[1]+".tmp")
		
		except:
			print '''Please check that misa.pl and misa.ini are located in the same 
			folder as this script'''
 
		try:
			finput = open(sys.argv[1]+".tmp", "rU")
			minput = open(sys.argv[1]+".tmp.misa", "rU")
		#	fastas = fastareader(finput)
		#	ssrs = misareader(minput)
		#	outputcreator(fastas,ssrs)

		except: 
			print "Could not open file. Please check location and read permissions"
		
		fastas = fastareader(finput)
		ssrs = misareader(minput)
		outputcreator(fastas,ssrs)
		garbagecollection(sys.argv[1])

## CREATE SORTING-FUNCTION FOR DICT ## 

def hashkeygen(record):
	description = record.description
	return description

## GET SEQUENCES OUT OF FASTA-FILE ##

def fastareader(finput):
	records = SeqIO.to_dict(SeqIO.parse(finput, "fasta"),key_function=hashkeygen)
	finput.close()
	return records

## READ THE SSRS OUT OF MISA-OUTPUT ##

def misareader(minput):
	ssrs = minput.readlines()
	ssrs.pop(0)
	ssroutput = []
	for i in ssrs:
		single_repeat = i.rstrip()
		single_repeat = single_repeat.replace("_", " ")
		single_repeat = single_repeat.split("\t")
		ssroutput.append(single_repeat)
	return ssroutput

## WRITE OUTPUT-FILE ##

def outputcreator(fastas, ssrs):
	outputliste = []
	for i in ssrs:
		outstring = i[0]+";"+i[1]+";"+i[2]+";"+i[3]+";"+i[4]+";"+i[5]+";"+i[6]+";"+str(fastas[i[0]].seq)+"\n"

		outputliste.append(outstring)
	output_handle = open(sys.argv[1]+".ssrs.csv","w")
	output_handle.write("ID;SSR nr.;SSR type;SSR;size;start;end;sequence\n")
	for i in outputliste:
		output_handle.write(i)
	output_handle.close()
	print "Finished. Output can be found at "+sys.argv[1]+".ssrs.csv"

## GET RID OF THOSE NASTY-FORMATTING-THINGS IN FASTA-FILES ##

def descriptorunifier(infile):
	inputfile = open(infile, "r")
	originaldata = inputfile.readlines()
	tempdata = []
	for i in originaldata:
		if i.find(">") == 0:
			i = i.replace("_", " ")
		tempdata.append(i)
	outfile = open(infile+".tmp","w")
	for i in tempdata:
		outfile.write(i)

## (RE-)MOVE TMP-FILES ## 

def garbagecollection(infile):
	os.system("rm "+infile+".tmp")
	os.system("mv "+infile+".tmp.misa "+infile+".misa")
	os.system("mv "+infile+".tmp.statistics "+infile+".statistics")

main()
