#! /usr/bin/env python

#import libraries for use
import re
import sys
import dendropy
import simple_matrix
import gene
import printer

#create list of trees for dendropy tree objects
treelist=dendropy.TreeList()
#list of genes, a class containing start and finish points for each gene
#allows for inter- and intragene comparisons
genelist=[]
#variable used to assign gene start and finish points
placeholder=0
#current gene name, used for transitioning between genes
CG=""
#reassign argv to not include 1st item (doesn't work for regex)
sys.argv=sys.argv[1:]
geneCount=0
for arg in sys.argv:
	#search file name for tree metadata
	#gene is group 1, run is group 2
	searchTerm="Alltaxa_(\w+)_r(\d)"
	regexobj=re.search(searchTerm,arg)
	#if the file doesn't match the current gene, make a new one and add it to genelist
	if regexobj.group(1) != CG:
		treelist.read_from_path(arg,'nexus')
		genelist.append(gene.gene(regexobj.group(1),placeholder,len(treelist)-1))
		geneCount+=1
	#if the next file is still the same gene, then extend finish on the gene
	else:
		genelist[-1].setFinish(len(treelist)-1)
		treelist.read_from_path(arg,'nexus')
		
	#assign metadata to the trees added from the current file specified in argv
	x=0
	for tree in treelist[placeholder:]:
		tree.name=regexobj.group(1)+"tree_"+str(x)
		x+=1
		tree.gene=regexobj.group(1)
		tree.geneNum=geneCount
		tree.run=regexobj.group(2)
		tree.freqCount=1
	placeholder=len(treelist)

#create distance matrix
Neo=simple_matrix.distMatrix(len(treelist))
#do the intragene tree comparisons first
print "num trees pre-pruning:",len(treelist)
for x in range(len(genelist)):
	#begin at gene.start, end at finish
	rc=genelist[x].start
	#the +1 is necessary because it stops 1 too early otherwise
	#TODO change finish indexing
	while rc < genelist[x].finish+1:
		#making cc larger than rc removes the self-comparison in the matrix
		cc=rc+1
		while cc < genelist[x].finish+1:
			rf=dendropy.treecalc.symmetric_difference(treelist[rc],treelist[cc])
			#if the distance isn't 0, then the value should be added to the matrix
			if rf!= 0:
				Neo.setVal(rc,cc,rf)
			else:
				#there's one more identical topology to tree at rc, increase freq
				#then remove tree from matrix and list
				treelist[rc].freqCount+=1
				Neo.removeIndex(cc)
				treelist.pop(cc)
				#also modify gene starts and finishes
				genelist[x].finish-=1
				for gen in genelist[x+1:]:
					gen.start-=1
					gen.finish-=1
			cc+=1
		#calculate tree frequency only after it's been compared to all other trees in gene
		treelist[rc].frequency=treelist[rc].freqCount/(float(genelist[x].totalsize))
		rc+=1
#debug prints to determine whether data fields are being assigned properly
"""for tree in treelist:
	print tree.name
	print tree.gene
	print tree.freqCount
	print tree.frequency"""
print "num trees post-pruning:",len(treelist)

#intergene comparison; something here is screwed up
for x in range(len(genelist)):
	#begin at gene.start, end at finish
	rc=genelist[x].finish+1
	
	#the +1 is necessary because it stops 1 too early otherwise

	
	#TODO change finish indexing
	while rc < len(treelist):
		#making cc larger than rc removes the self-comparison in the matrix
		cc=rc+1
		while cc < len(treelist):
			#print "rc:",rc
			#print "name:",treelist[rc].name
			#print "cc:",cc
			#print "name:",treelist[cc].name
			rf=dendropy.treecalc.symmetric_difference(treelist[rc],treelist[cc])
			#print "rf:",rf
			Neo.setVal(rc,cc,rf)
			#print "matrix value:",Neo.rowlist[rc][cc]
			cc+=1
		rc+=1
		
myPrint = printer.printManager(treelist,genelist,Neo)

for x in range(len(genelist)-2):
	for y in range(x+1,len(genelist)-1):
		for z in range(y+1,len(genelist)):
			myPrint.JhiveOutput([x,y,z],.02)


myPrint.JhiveOutput([1,2,3],.02)
