#! /usr/bin/env python

import re
import sys
import dendropy
import simple_matrix
import gene

class printManager:
	def __init__(self,treelist,genelist,matrix):
		self.treelist=treelist
		self.genelist=genelist
		self.Neo=matrix

	def HivePlotterOutput(self,genewrite,threshold):
		#genewrite; maybe make loop to iterate different combos?
		genedictionary={0:'a',1:'b',2:'c',3:'d',4:'e'}
		#HTML node output
		nodename="HTML_node_"+str(genewrite[0])+str(genewrite[1])+str(genewrite[2])+".txt"
		output=open(nodename,'w')
		output.write("Node\tFrequency\tGene\n")
		for x in range(len(genewrite)):
			for tree in self.treelist[self.genelist[genewrite[x]].start:self.genelist[genewrite[x]].finish]:
				if tree.frequency > threshold:
					output.write(tree.name+"\t"+str(tree.frequency) +"\t"+tree.gene+"\n")

		output.close()

		#HTML edge output
		edgename="HTML_edge_"+str(genewrite[0])+str(genewrite[1])+str(genewrite[2])+".txt"
		output=open(edgename,'w')
		output.write("Source\ttarget\tDistance\n")

		#expand comments to more fully explain comparisons, could get weird for more than 3
		#go through each gene in genewrite, set the y coord as that gene, and x as the next
		for x in range(len(genewrite)-1):
			for y in range(self.genelist[genewrite[x]].start,self.genelist[genewrite[x]].finish):
				for z in range(self.genelist[genewrite[x+1]].start,self.genelist[genewrite[x+1]].finish):
					#print "z:",z
					if self.treelist[y].frequency > .03 and self.treelist[z].frequency > .03:
						output.write(self.treelist[y].name+"\t"+self.treelist[z].name+"\t"+str(self.Neo.rowlist[y][z])+"\n")

		for x in range(self.genelist[genewrite[-1]].start,self.genelist[genewrite[-1]].finish):
			for y in range(self.genelist[genewrite[0]].start,self.genelist[genewrite[0]].finish):
				if self.treelist[x].frequency > .03 and self.treelist[y].frequency > .03:
					output.write(self.treelist[x].name+"\t"+self.treelist[y].name+"\t"+str(self.Neo.rowlist[y][z])+"\n")
		output.close()
		
	def JhiveOutput (self,genewrite,threshold):
		filename="Jhive_"+str(genewrite[0])+"_"+str(genewrite[1])+"_"+str(genewrite[2])+".dot"
		output=open(filename,"w")
		#write trees first, per dot format
		for x in range(len(genewrite)):
			for tree in self.treelist[self.genelist[genewrite[x]].start:self.genelist[genewrite[x]].finish]:
				if tree.frequency > threshold:
					# output looks like
					#ATP6tree_0[ name=ATP6tree_0 gene=ATP6 geneNumber=1 frequency=0.18]
					output.write(tree.name+"[ name="+tree.name +" gene="+tree.gene +" geneNumber="+str(x) + " frequency=" + str(tree.frequency)+"]\n")
		
		#write edges next, in same file		
		for x in range(len(genewrite) - 1):
			for y in range(self.genelist[genewrite[x]].start,self.genelist[genewrite[x]].finish):
				for z in range(self.genelist[genewrite[x+1]].start,self.genelist[genewrite[x+1]].finish):
					if self.treelist[y].frequency > threshold and self.treelist[z].frequency > threshold:
						output.write(self.treelist[y].name + "--" + self.treelist[z].name + " [rf="+str(self.Neo.rowlist[y][z])+"]\n")
						
		for x in range(self.genelist[genewrite[-1]].start,self.genelist[genewrite[-1]].finish):
			for y in range(self.genelist[genewrite[0]].start,self.genelist[genewrite[0]].finish):
				if self.treelist[x].frequency > threshold and self.treelist[y].frequency > threshold:
					output.write(self.treelist[x].name + "--" + self.treelist[y].name + " [rf="+str(self.Neo.rowlist[x][y])+"]\n")
		output.close()
