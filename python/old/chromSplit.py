#!/usr/bin/env python

#Import two bed files
#	1.summary BED(chr,start,stop,fold_change)
#	2.gene BED(chr,start,stop,id)
#Set up export for each chromosome
#For each line in each chromosome block of summary BED,
#	 find interval in gene bed such that:
#		summary.start>gene.start-1000
#		AND
#		summary.stop<gene.start+1000
#	 write to export (id,value)

import re, os, shutil, time, sys
from math import *
from string import *
from optparse import OptionParser
import operator

import BED_mod
import GenomeData
import lineCount

"""
Count number of lines in gene file
For each test
        bisect gene file
        if intersect: write
        elif greater: bisect above
        elif lesser: bisect below
"""


def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")



	
	(opt, args) = parser.parse_args(argv)
	
	input = open(opt.input,'r')
	tn = ""
	tline = input.readline();
	stline = tline.split();
        bed_length = len(stline)
        #sys.exit();
        types = "";
        
        if bed_length == 4: types = "BED_GRAPH"
                        
        elif bed_length== 6: types ="BED6"
        
	if type(stline[3]) == type("x"):
		tn = "name"
        
	
	testBED = BED_mod.BED(species="mm9",file=opt.input,bed_type=types,n=tn);
	outfiles = dict()
	for chrom in testBED.keys():
		file = opt.output + "_" + chrom 
		outfiles[chrom] = file

	print outfiles;
	#sys.exit(); 
	
	for chrom in testBED.keys():
                #count = len(geneBED.bed_vals[chrom]);
                #sys.exit();
		#print "Test chromosome: " + chrom;
		#print len(geneBED.bed_vals[chrom]);
		#print geneBED.bed_vals[chrom][1].start
		out = open(outfiles[chrom],'w');
		for val in testBED.bed_vals[chrom]:
			if types == "BED6":
				o = chrom + "\t" + str(val.start) + "\t" + str(val.end) + "\t" + str(val.name) + "\t" + str(val.score) + "\t" + str(val.strand) + "\n"
			elif types == "BED_GRAPH":
				o = chrom + "\t" + str(val.start) + "\t" + str(val.end) + "\t" + str(val.value) + "\n"
			out.write(o)
		out.close();


def conditions(ts,te,gs,ge,o,s):
	a = 0;
	if ts<=gs and te>=ge:
		a = 1;	
	elif ts>=gs and te<=ge: 
		if s:	a = 1;
		else:	a = 2;
	elif ts<=gs and te>=gs+o:
		if s:	a = 1;
		else:	a = 3;
	elif ts<=ge+o and te>=ge: 
		if s:	a = 1;
		else:	a = 4;
	return a; 	
	
			
	
	
	
def conditions2(ts,te,gs,ge,o,f):
	a = 0
	if ts<=gs and te>=ge:		
		a = 1;
	elif ts>=gs and te<=ge:
		a=2;
	elif ts<=gs and te>=gs+o:
		a=3;
	elif ts<=ge+o and te>=ge:
		a=4;
	return a;
			

if __name__ == "__main__":
	main(sys.argv)
