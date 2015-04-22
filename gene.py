#! /usr/bin/env python

class gene:
	"""
	class representing the start and ending positions of trees belonging to a particular gene
	dynamic start and finish points allow for trees to be added or removed
	"""
	def __init__(self,name,start,finish):
		self.name=name
		self.start=start
		self.finish=finish
		self.totalsize=(finish+1)-start
	def setStart(num):
		"""
		sets the start index of gene object to a specified number
		"""
		self.start=num
		self.totalsize=(finish+1)-start
	def setFinish(num):
		"""
		sets the finish index of gene object to a specified number
		"""
		self.finish=num
		self.totalsize=(finish+1)-start
	def shiftStart(shift):
		"""
		moves the start index backward or forward by adding the argument to the original index
		"""
		self.start=self.start + shift
		self.totalsize=(finish+1)-start
	def shiftFinish(shift):
		"""
		moves the finish index backward or forward by adding the argument to the original index
		"""
		self.finish=self.finish + shift
		self.totalsize=(finish+1)-start
