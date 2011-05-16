#!/usr/bin/python

from Bio.Blast import NCBIXML
import sys

# TRY TO OPEN THE INPUT FILE #

def input():
	try:
		infile = open(sys.argv[1])
		return infile
	except: 
		print "This script takes a standard NCBI-Blast XML-file and converts it"
		print "into a tab-seperated outfile which can be opened using standard"
		print "office-tools like Open Office Calculator or MS Excel"
		print ""
		print "Just start the script giving a single XML-File"
		print "i.e. \"./xmlextractor.py blastout.xml\""


# GET DATA OUT OF XML-FILE #
# DATA WILL BE STORED IN NESTED LIST [[identifier, foo, bar], [identifier, foo, bar],]]
def parser(handler):
	blast_records = NCBIXML.parse(handler)
	outlist = []
	outlist.append("Read\tAlignment Title\tScore\tBits\teValue\tIdentity+\tStart\tStop\n")
	for i in blast_records:
		read = i.query
		for alignment in i.alignments:
			aligntitle = alignment.title
			for hsp in alignment.hsps:
				score = hsp.score
				bits = hsp.bits
				expect = hsp.expect
				identity = hsp.identities
				start = hsp.query_start
				query = len(hsp.query*3)
				if len(sys.argv) < 3:
					outlist.append(read+"\t"+aligntitle+"\t"+str(score)+"\t"+str(bits)+"\t"+str(expect)+"\t"+str(identity)+"\t"+str(start)+"\t"+str(start+query-1)+"\n")
				elif expect < float(sys.argv[2]):
					outlist.append(read+"\t"+aligntitle+"\t"+str(score)+"\t"+str(bits)+"\t"+str(expect)+"\t"+str(identity)+"\t"+str(start)+"\t"+str(start+query-1)+"\n")
	return outlist		

def writer(input):
	outfile = open(sys.argv[1]+".tab","w")
	for i in input:
		outfile.write(i)

handler = input()
if handler != None:
	outlist = parser(handler)
	writer(outlist)
