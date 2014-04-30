#!/usr/bin/env python

import sys
import os
import math
import copy
import numpy as np
import matplotlib.pyplot as plt


class Interval:
	def __init__(self):
		self.chr = ""
		self.start = -1
		self.end = -1
		self.sig = 0.00
		return

class Bed:
	def __init__(self):
		self.intervals = []
		return

	def sortByChromosomeAndStartAndEnd(self):
		self.intervals.sort(key = lambda Interval:Interval.end)
		self.intervals.sort(key = lambda Interval:Interval.start)
		self.intervals.sort(key = lambda Interval:Interval.chr)
		#chrList = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
		#intervalList = []
		#for eachChr in chrList:
			#tempIntervals = []
			#for currInterval in self.intervals:
				#if currInterval.chr == eachChr:
					#copyInterval = Interval()
					#copyInterval.chr = currInterval.chr
					#copyInterval.start = currInterval.start
					#copyInterval.end = currInterval.end
					#copyInterval.sig = currInterval.sig
					#tempIntervals.append(copyInterval)
			#tempIntervals.sort(key = lambda Interval:Interval.end)
			#tempIntervals.sort(key = lambda Interval:Interval.start)
			#for eachInterval in tempIntervals:
				#copyInterval = Interval()
				#copyInterval.chr = eachInterval.chr
				#copyInterval.start = eachInterval.start
				#copyInterval.end = eachInterval.end
				#copyInterval.sig = eachInterval.sig
				
				#intervalList.append(copyInterval)
		#self.intervals = []
		#self.intervals = intervalList	
		#for eachInterval in self.intervals:
			#print eachInterval.chr, eachInterval.start, eachInterval.end
		return

class Region:
	def __init__(self):
		self.chr = ""
		self.start1 = -1
		self.end1 = -1
		self.start2 = -1
		self.end2 = -1
		return

def readBedPeaks(filename):
	ip = open(filename, "r")
	currBed = Bed()
	avgPeak = 0
	for line in ip:
		if ("#" not in line) and (len(line)) >= 3:
			fields = line.split()
			currInterval = Interval()
			currInterval.chr = fields[0]
			currInterval.start = int(fields[1])
			currInterval.end = int(fields[2])
			avgPeak += currInterval.end - currInterval.start
			currBed.intervals.append(currInterval)
	avgPeak = float(avgPeak)/len(currBed.intervals)
	print "averageWidth is", avgPeak
	ip.close()
	return currBed

def greaterThan(interval1, interval2):
	if interval2.chr > interval1.chr:
		return True
	elif interval1.chr == interval2.chr and interval2.start > interval1.end:
		return True
	else:
		return False


def skipToCurrInterval(ipSig, currInterval):
	line = ipSig.readline()
	if line.split()[0] > currInterval.chr:
		ipSig.seek(0)
		line = ipSig.readline()
	
	intersecting = False
	while not(intersecting):
		#print currInterval.chr, currInterval.start, currInterval.end, line
		if ("#" not in line) and (len(line)) >= 3:
			fields = line.split()
			testInterval = Interval()
			testInterval.chr = fields[0]
			testInterval.end = int(fields[2])
			if testInterval.chr == currInterval.chr and testInterval.end >= currInterval.start: 
				testInterval.start = int(fields[1])
				testInterval.sig = float(fields[3])
				intersecting = True
				break
		line = ipSig.readline()
	return testInterval

def readBedGraphForPeak(ipSig, currInterval):
	currBed = Bed()
	line = ipSig.readline()
	intersectingInterval = skipToCurrInterval(ipSig, currInterval)
	#print currInterval.chr, currInterval.start, currInterval.end, intersectingInterval.chr, intersectingInterval.start, intersectingInterval.end
	currBed.intervals.append(intersectingInterval)
	intersecting = True
	while intersecting:
		line = ipSig.readline()
		if ("#" not in line) and (len(line)) >= 4:
			fields = line.split()
			intersectingInterval = Interval()
			intersectingInterval.chr = fields[0]
			intersectingInterval.start = int(fields[1])
			intersectingInterval.end = int(fields[2])
			intersectingInterval.sig = float(fields[3])
			if (intersectingInterval.start > currInterval.end):
				#print currInterval.chr, currInterval.start, currInterval.end, intersectingInterval.chr, intersectingInterval.start, intersectingInterval.end
				intersecting = False
				break
			currBed.intervals.append(intersectingInterval) 
	return currBed

def smoothSignal(bedSig, width):
	idx = width/2
	length = len(bedSig.intervals)
	while width/2 <= idx <= (length - (width/2) - 1):
		currSig = 0.
		idx2 = idx - width/2
		while idx2 <= idx + width/2:
			currSig += bedSig.intervals[idx2].sig
			idx2 += 1
		currSig /= width
		bedSig.intervals[idx].sig = currSig
		idx += 1
	return

def findPeaks(bedSig, peaks, lowWidth, width):
	
	#Checking for peak in first few bins
	length = len(bedSig.intervals)
	if length < width/2:
		return
	signals = []
	for eachInterval in bedSig.intervals:
		signals.append(eachInterval.sig)
	maximumSig = max(signals)
	#print maximumSig
	for idx in range(0, width/2):
		if math.fabs(maximumSig - bedSig.intervals[idx].sig) < 0.001:
			currInterval = Interval()
			currInterval.chr = bedSig.intervals[idx].chr
			currInterval.start = bedSig.intervals[idx].start
			currInterval.end = bedSig.intervals[idx].end
			currInterval.sig = bedSig.intervals[idx].sig
			
			if len(peaks.intervals) > 1 and  bedSig.intervals[idx].chr == peaks.intervals[-1].chr and bedSig.intervals[idx].start < peaks.intervals[-1].end + lowWidth:
				continue
			peaks.intervals.append(currInterval)
			#print currInterval.chr, currInterval.start, currInterval.end
			break

	#Checking for peak in middle bins
	idx = width/2
	while width/2 <= idx <= (length - (width/2) - 1):
		currSig = bedSig.intervals[idx].sig
		idx2 = idx - width/2
		while idx2 <= idx + width/2:
			if currSig < bedSig.intervals[idx2].sig:
				break
			idx2 += 1
		if idx2 > (idx + width/2) and currSig >= maximumSig/3:
			if len(peaks.intervals) > 1 and  bedSig.intervals[idx].chr == peaks.intervals[-1].chr and bedSig.intervals[idx].start < peaks.intervals[-1].end + lowWidth:
				idx += 1
				continue
			currInterval = Interval()
			currInterval.chr = bedSig.intervals[idx].chr
			currInterval.start = bedSig.intervals[idx].start
			currInterval.end = bedSig.intervals[idx].end
			currInterval.sig = bedSig.intervals[idx].sig
			peaks.intervals.append(currInterval)
			#print peaks.intervals[-1].chr, peaks.intervals[-1].start, peaks.intervals[-1].end, peaks.intervals[-1].sig
		idx += 1

	#Checking for peak in last few bins
	for idx in range((length - (width/2)), length):
		if math.fabs(maximumSig - bedSig.intervals[idx].sig) < 0.001:
			currInterval = Interval()
			currInterval.chr = bedSig.intervals[idx].chr
			currInterval.start = bedSig.intervals[idx].start
			currInterval.end = bedSig.intervals[idx].end
			currInterval.sig = bedSig.intervals[idx].sig
			if len(peaks.intervals) > 1 and  bedSig.intervals[idx].chr == peaks.intervals[-1].chr and bedSig.intervals[idx].start < peaks.intervals[-1].end + lowWidth:
				continue
			peaks.intervals.append(currInterval)
			#print currInterval.chr, currInterval.start, currInterval.end
			break
	return


def checkForPeaksWithinWidth(regionList, histonePeaks, minWidth, maxWidth):
	length = len(histonePeaks.intervals)
	for idx1 in range(0, length):
		center1 = (histonePeaks.intervals[idx1].end + histonePeaks.intervals[idx1].start)/2
		for idx2 in range(idx1 + 1, length):
			center2 = (histonePeaks.intervals[idx2].end + histonePeaks.intervals[idx2].start)/2
			if (histonePeaks.intervals[idx1].chr == histonePeaks.intervals[idx2].chr) and ((center2 - center1) > minWidth) and ((center2 - center1) <= maxWidth):
				currRegion = Region()
				currRegion.chr = histonePeaks.intervals[idx1].chr 
				currRegion.start1 = histonePeaks.intervals[idx1].start
				currRegion.end1 = histonePeaks.intervals[idx1].end
				currRegion.start2 = histonePeaks.intervals[idx2].start
				currRegion.end2 = histonePeaks.intervals[idx2].end
				regionList.append(currRegion)
				#print currRegion.chr, currRegion.start1, currRegion.end1, currRegion.start2, currRegion.end2
			elif (histonePeaks.intervals[idx1].chr != histonePeaks.intervals[idx2].chr):
				break
			elif ((center2 - center1) > maxWidth):
				break
			idx2 += 1
		idx1 += 1
	return


def mergeOverlappingRegions(regions):
	length = len(regions)
	idx1 = 0
	newRegions = []
	while idx1 < length:
		currRegion = Region()
		currRegion.chr = regions[idx1].chr
		currRegion.start1 = regions[idx1].start1
		currRegion.end1 = regions[idx1].end1
		currRegion.start2 = regions[idx1].start2
		currRegion.end2 = regions[idx1].end2
		idx2 = idx1 + 1
		while idx2 < length:
			if currRegion.chr == regions[idx2].chr and currRegion.start1 >= regions[idx2].start1 and currRegion.start1 <= regions[idx2].end2 and currRegion.end2 >= regions[idx2].end2:
				#print currRegion.chr, currRegion.start1, currRegion.end2, regions[idx2].start1, regions[idx2].end2
				currRegion.start1 = regions[idx2].start1
				currRegion.end1 = regions[idx2].end1
				idx2 += 1
			elif currRegion.chr == regions[idx2].chr and currRegion.start1 <= regions[idx2].start1 and currRegion.end2 >= regions[idx2].end2:
				#print currRegion.chr, currRegion.start1, currRegion.end2, regions[idx2].start1, regions[idx2].end2
				idx2 += 1
			elif currRegion.chr == regions[idx2].chr and currRegion.start1 <= regions[idx2].start1 and currRegion.end2 >= regions[idx2].start1:
				currRegion.start2 = regions[idx2].start2
				currRegion.end1 = regions[idx2].end2
				idx2 += 1
			else:
				break		
		idx1 = idx2 + 1
		newRegions.append(currRegion)
	#print len(newRegions)
	idx = len(regions)
	while idx > 0:
		regions.pop()
		idx -= 1
	for currRegion in newRegions:
		regions.append(currRegion)
	#regions = copy.deepcopy(newRegions)
	#print len(newRegions), len(regions)
	pairedWidth = 0
	for currRegion in regions:
		pairedWidth += currRegion.end2 - currRegion.start1
	pairedWidth = float(pairedWidth)/len(regions)
	print "pairedWidth =", pairedWidth
	return

def sortRegionsByChromosomeAndStartAndEnd(regions):
	regions.sort(key = lambda Regions:Regions.end2)
	regions.sort(key = lambda Regions:Regions.start1)
	regions.sort(key = lambda Regions:Regions.chr)
	return

def checkForPeaksWithinRegion(chromPeaks, doubleConditionRegions, regionList):
	length = len(chromPeaks.intervals)
	idx = 0
	for currRegion in regionList:
		#print currRegion.chr, currRegion.start1, currRegion.end2
		#print chromPeaks.intervals[idx].chr, chromPeaks.intervals[idx].start, chromPeaks.intervals[idx].end
		while idx < len(chromPeaks.intervals):
			if currRegion.chr == chromPeaks.intervals[idx].chr and (currRegion.start1 - 10000) <= chromPeaks.intervals[idx].end and (currRegion.end2 + 10000) >= chromPeaks.intervals[idx].start:
				if currRegion not in doubleConditionRegions:
					doubleConditionRegions.append(currRegion)
					#print idx, currRegion.chr, currRegion.start1, currRegion.end2, chromPeaks.intervals[idx].start, chromPeaks.intervals[idx].end
				break
			if currRegion.chr < chromPeaks.intervals[idx].chr:
				break
			if currRegion.end2 < chromPeaks.intervals[idx].start:
				break
			idx += 1
	return


if (len(sys.argv) != 3):
	sys.stderr.write("Usage: " + sys.argv[0] + " <histoneInputFileList> <outputBedFile>\n")
	sys.exit()

lowWidth = 400 #min of 0.2kb width between 2 peaks
highWidth = 2000 #max of 2kb width between 2 peaks
ipHistone = open(sys.argv[1], "r")
regionList = []
finalHistonePeaks = Bed()
for line in ipHistone:
	histonePeaks = readBedPeaks(line.rstrip().split("\t")[0])
	histonePeaks.sortByChromosomeAndStartAndEnd()
	ipSig = open(line.rstrip().split("\t")[1], "r")
	for currInterval in histonePeaks.intervals:
		#print currInterval.chr, currInterval.start, currInterval.end
		histoneSig = readBedGraphForPeak(ipSig, currInterval)
		sigRaw = []
		posn = []
		for currInterval in histoneSig.intervals:
			sigRaw.append(currInterval.sig)
			posn.append(currInterval.start)
		#print len(histoneSig.intervals)
		smoothSignal(histoneSig, 5)
		sigSmooth = []
		for currInterval in histoneSig.intervals:
			sigSmooth.append(currInterval.sig)
		findPeaks(histoneSig, finalHistonePeaks, lowWidth, 5)
		#if currInterval.chr == "chr10":
		#print "Got back"
		#plt.plot(posn, sigRaw, posn, sigSmooth)
		#plt.show()
		#sys.exit()
	print len(finalHistonePeaks.intervals)
	checkForPeaksWithinWidth(regionList, finalHistonePeaks, lowWidth, highWidth)
	print len(regionList)
	sortRegionsByChromosomeAndStartAndEnd(regionList)
	mergeOverlappingRegions(regionList)
	print len(regionList)
ipHistone.close()
#sys.exit()

sortRegionsByChromosomeAndStartAndEnd(regionList)
#doubleConditionRegions = []
#ipChrom = open(sys.argv[2], "r")
#for line in ipChrom:
	##chromPeaks = readBedPeaks(line.rstrip())
	#chromPeaks.sortByChromosomeAndStartAndEnd()
	#checkForPeaksWithinRegion(chromPeaks, doubleConditionRegions, regionList)
	#print len(doubleConditionRegions)
	#ipSig = open(line.rstrip())
	#for currRegion in regionList:
		#currChromSig = readWiggleCurrInterval(ipSig, currRegion.chr,  currRegion.start1, currRegion.end2)


#ipChrom.close()

op = open(sys.argv[2], "w")
for currRegion in regionList:
	op.write(currRegion.chr + "\t" + str(currRegion.start1) + "\t" + str(currRegion.end2) + "\n")
op.close()
