#! /usr/bin/env python
import sys
import numpy as np
import re
import copy

class bipartition:
	def __init__(self,name,code,freq):
		self.name,self.code,self.freq=name,code,freq
		
class community:
	def __init__(self):
		self.bipartitionDictionary={}

print "script started"
#command line arguments
# first file is community file
#second file is covariance matrix
filenameOne=sys.argv[1]
filenameTwo=sys.argv[2]

infile=open(filenameOne,'r')

#taxon dictionary used to retrieve taxon name from index number
taxDict={}
#bipartition dictionary
bipartDict={}
#community dictionary
communityDict={}
#list to store the number value of "floating" bipartitions
floaters=[]

isTaxa=False
#regular expression strings to pull out necessary info
taxSearch="(t\d+s\d+)\;* , (\d+)"
bipartSearch="bipartition (\d+) : (\d+), appear times: (\d+)"
communitySearch="Community (\d+) includes nodes: ((\d+,)+)"

#skips all the crap in the file prior to the taxa listed; sets up for next loop and taxon info
while isTaxa==False:
	line=infile.readline()
	print "taxaSearch:",line
	regex=re.search(taxSearch,line)
	if regex == None:
		pass
	else:
		isTaxa=True
		taxDict[int(regex.group(2))]=regex.group(1)
		print "notTaxa done"

#now go through the taxa until bipartitions come up
while isTaxa:
	line=infile.readline()

	regex=re.search(taxSearch,line)
	
	if regex != None:
		taxDict[int(regex.group(2))]=regex.group(1)
	else:
		isTaxa=False
		regex=re.search(bipartSearch,line)
		bipartDict[int(regex.group(1))]=bipartition.bipartition(regex.group(1),regex.group(2),regex.group(3))
		print "isTaxa done"

#create bipartition objects and add them to list
isBipart=True
while isBipart:
	line=infile.readline()
	regex=re.search(bipartSearch,line)
	
	if regex !=None:
		bipartDict[int(regex.group(1))]=bipartition.bipartition(regex.group(1),regex.group(2),regex.group(3))
	else:
		isBipart=False
		print "isBipart done"

#skip a bunch of crap in the file; maybe save some later if need be
isCommunity=False
while isCommunity==False:
	line=infile.readline()
	regex=re.search(communitySearch,line)
	if regex == None:
		pass
	else:
		isCommunity=True
		numBiparts=regex.group(2)
		numBiparts=numBiparts[:-1]
		bipartList=numBiparts.split(",")
		communityDict[int(regex.group(1))]=community.community()
		for number in bipartList:
			#take out the +1 in newer files since indexing issues went away
			communityDict[int(regex.group(1))].bipartitionDictionary[int(number)+1]=bipartDict[int(number)+1]
		print "notCommunity done"

#go through various communities and assign different bipartitions to them
while isCommunity:
	line=infile.readline()
	regex=re.search(communitySearch,line)
	if regex != None:
		numBiparts=regex.group(2)
		numBiparts=numBiparts[:-1]
		bipartList=numBiparts.split(",")
		communityDict[int(regex.group(1))]=community.community()
		for number in bipartList:
			#take out the +1 in newer files since indexing issues went away
			communityDict[int(regex.group(1))].bipartitionDictionary[int(number)+1]=bipartDict[int(number)+1]
	else:
		isCommunity = False
		print "isCommunity done"

#string containing floating bipartitions
line=infile.readline()
#get rid of annoying spaces
line=line.replace(" ","")	
#remove new line character as well as floating comma
line=line[0:-2]
#get list of bipartitions from line and add to floaters list; a list of strings, change to ints?
floaterIndices=line.split(",")
floaters.extend(floaterIndices)

#last few lines of file don't matter; close it!		
infile.close()

#go through communities to find biggest
#set initial communities to nothing
commOneSize=0
commOneIndex=-1
commTwoSize=0
commTwoIndex=-1
bigCommTwo=	-1
#all i's are +1 due to stupid indexing
for i in range(len(communityDict)):

	if len(communityDict[i+1].bipartitionDictionary) >= commOneSize:
		#set second largest from current largest
		commTwoSize=commOneSize
		commTwoIndex=commOneIndex
		#reassign largest community value
		commOneSize=len(communityDict[i+1].bipartitionDictionary)
		commOneIndex=i+1
		


#read in pairwise-covariance matrix
infile=open(filenameTwo,'r')
infile.readline()
index=0
redPill=[]

for line in infile:
	#get rid of newline
	line=line.strip()
	line=line.replace(" ","")
	#print "line:", line
	#tab delimited
	vals=line.split("\t")
	#remove leading row label
	vals=vals[1:]
	#print "val:", vals
	redPill.append(vals)	
	#print "list o' lists:",myMatrix
	#increment index	
	index+=1

#myMatrix is non-symmetrical, need to fill in a complete matrix

bluePill=np.matrix(np.empty((len(redPill),len(redPill))))

for i in range(len(bluePill)):
	for j in range(i,len(bluePill)):
			#assign values across the diagonal!
			#print i,j
			#print redPill[j]
			bluePill[i,j]=redPill[j][i]
			bluePill[j,i]=redPill[j][i]

#make output file

#all the i+1 stuff can be changed to i in files w/o indexing error of bipartitions (i.e. actually start w/ 0)
outfile=open("Jhive_bipart.dot",'w')
for i in range(len(bipartDict)):
	#if the community is found in a community, that will effect the variable used to sort by axis
	if i in communityDict[commOneIndex].bipartitionDictionary.keys():
		outfile.write(bipartDict[i+1].name + "[ name=" + bipartDict[i+1].name +" frequency="+str(bipartDict[i+1].freq)+ " axis=1]\n")
	elif i in communityDict[commTwoIndex].bipartitionDictionary.keys():
		outfile.write(bipartDict[i+1].name + "[ name=" + bipartDict[i+1].name +" frequency="+str(bipartDict[i+1].freq)+ " axis=2]\n")
	elif str(i) in floaters:
		outfile.write(bipartDict[i+1].name + "[ name=" + bipartDict[i+1].name +" frequency="+str(bipartDict[i+1].freq)+ " axis=3]\n")
	else:
		outfile.write(bipartDict[i+1].name + "[ name=" + bipartDict[i+1].name +" frequency="+str(bipartDict[i+1].freq)+ " axis=0]\n")
#now that nodes and edges are in memory, start writing stuff!
#write nodes to file, changing axis assignment variable based on whether node belongs to largest communities
covarMax=0
#writes edge info; find min and max covariance values
for x in range(len(bipartDict)-1):
	for y in range(x,len(bipartDict)-1):
		for z in range(x+1,len(bipartDict)-1):
			#the +1 is necessary since communities aren't 0 indexed
			if bluePill[y,z] != 0:
				outfile.write(bipartDict[y+1].name + "--" + bipartDict[z+1].name + "[covar="+str(bluePill[y,z])+ "]\n")
				#check matrix value, find maximum absolute covariance value
				if abs(bluePill[y,z]) > covarMax:
					covarMax=abs(bluePill[y,z])
outfile.close()

#formatting file stuff;colors edges according to an 11-category divergent heatmap
#default color is 0, reds are positive, and blues are negative

#need max and min covariance values
outfile=open("Jhive_bipart_format.txt",'w')
covarInc=covarMax/5
negColors=["(209,229,240)","(146,197,222)","(67,147,195)","(33,102,172)","(5,48,97)"]
posColors=["(253,219,199)","(244,165,130)","(214,96,77)","(178,24,43)","(103,0,31)"]

#write negative color values
for i in range(0,5):
	outfile.write("e(covar <="+str((i*-covarInc))+") [color="+negColors[i]+"]\n")
#write positive color values
for i in range(0,5):
	outfile.write("e(covar >="+str((i*covarInc))+") [color="+posColors[i]+"]\n")

outfile.close()
