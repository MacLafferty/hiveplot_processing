#! /usr/bin/env python

import sys
import re

#1st command-line argument 
filename=sys.argv[1]

infile=open(filename,'r')

#process various pre-taxa info
#
infile.readline()
infile.readline()
numTrees=infile.readline()

#pull out total number of bipartitions
numTotalBipartsString=infile.readline()
searchString="Unique bipartition number: (\d+)"
regexObj=re.search(searchString,numTotalBipartsString)
numTotalBiparts=int(regexObj.group(1))
infile.readline()


taxSearch="(t\d+s\d+) , (\d+)"
bipartSearch="bipartition (\d+) : (\d+), appear times: (\d+)"
taxDict={}
isTaxa=True
while isTaxa:
	line=infile.readline()
	regex=re.search(taxSearch,line)
	
	if regex != None:
		taxDict[regex.group(2)]=regex.group(1)
	else:
		isTaxa=False
		regex=re.search(bipartSearch,line)
