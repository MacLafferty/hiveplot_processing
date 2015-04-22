#! /usr/bin/env python
import dendropy


class distMatrix:
	def __init__(self,size):
		self.rowlist=[]
		for x in range(size):
			self.rowlist.append([])
			for y in range(size):
				self.rowlist[x].append(0)
			
	def setVal(self,x,y,val):
		self.rowlist[x][y]=val
		self.rowlist[y][x]=val
		
	def removeIndex(self,index):
		self.rowlist.pop(index)
		for x in range(len(self.rowlist)):
			self.rowlist[x].pop(index)

	def displayList(self):
		for x in range(len(self.rowlist)):
			print self.rowlist[x]
