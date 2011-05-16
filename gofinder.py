#!/usr/bin/python

import sys
from goatools import obo_parser

# a script that takes blast2go-output and finds names & namespaces to the go-ids in this data
# just pass a blast2go-outfile and the obo-database
print sys.argv

def goreader():
	infile = open(sys.argv[1],"r")
	gohash = obo_parser.GODag(sys.argv[2])
	outfile = open(sys.argv[1]+"_go.tab","w")

	inputdata = infile.readlines()
	splitdata = []
	outdata = []
	outfile.write("Read+\tGO-ID\tGO-Name\tGo-Namespace\n")

	for i in inputdata:
		splitdata.append(i.split())

	for j in splitdata:
		outfile.write(j[0]+"\t"+j[1]+"\t"+gohash[j[1]].name+"\t"+gohash[j[1]].namespace+"\n")


def main():
	try:
		x = open(sys.argv[1])
		y = open(sys.argv[2])
		goreader()
	except:
		
		print "This script converts standard BLAST2GO-Output into"
		print "somewhat more meaningful, tab-separated output"
		print "Besides the query and the resulting GO-ID the output"
		print "includes the name and the namespace of the GO-term"
		print ""
		print "To get started pass a BLAST2GO-file and a obo-file"
		print "which can be obtained from gene ontology"
		print ""
		print "i.e. \"./gofinder.py blast2go.annot geneontology.obo\""

main()
