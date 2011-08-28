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





def main(argv):

	parser = OptionParser()

	#parser.add_option("-s", "--species", action="store", type="string", dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-t", "--test", action="store", type="string", dest="test", metavar="<str>")
	parser.add_option("-g", "--gene", action="store", type="string", dest="gene", metavar="<str>")
	parser.add_option("-i", "--overlap", action="store", type="int", dest="overlap", metavar="<int>",default=0)
	parser.add_option("-o", "--outFile", action="store", type="string", dest="outFile", metavar="<file>")
	parser.add_option("-k", "--kind", action="store", type="string", dest="kind",metavar="<str>")
	parser.add_option("-n", "--name", action="store", type="string", dest="name",metavar="<str>")
	parser.add_option("-f", "--fraction", action="store", type="float", dest="fraction", metavar="<float>",default=0.0)

	
	(opt, args) = parser.parse_args(argv)
	
	overlap=opt.overlap;
	fraction=opt.fraction;
	#if opt.species in GenomeData.species_chroms.keys():
	#	print "Species: ", opt.species;
	gene = open(opt.gene,'r')
	test = open(opt.test,'r')
	gn = ""
	tn = ""
	gline,tline = gene.readline(),test.readline();
	sgline,stline = gline.split(),tline.split();

	if type(sgline[3]) == type("x"):
		gn = "name"
	if type(stline[3]) == type("x"):
		tn = "name"


	print "Accessing test BED...";
	testBED = BED_mod.BED(species="mm9",file=opt.test,bed_type="BED_GRAPH",n=tn);
	print "Accessing gene BED...";
	geneBED = BED_mod.BED(species="mm9",file=opt.gene,bed_type="BED_GRAPH",n=gn);
	outputfile = open(opt.outFile, 'w');
	total = 0;
	
	for chrom in testBED.keys():
		
		print "Test chromosome: " + chrom;
		if opt.kind == "in":
			for g in geneBED.bed_vals[chrom]:
				reads=0;
				switch=0;
				
				for t in testBED.bed_vals[chrom]:
					
					if conditions(t.start,t.end,g.start,g.end,overlap,1):
						reads += 1
						switch = 1;
					elif switch: 
						out = chrom + "\t" + str(g.start)+"\t"+str(g.end)+"\t"+str(g.value) + "\t" + str(reads) + "\n";
						outputfile.write(out);
						break;
	
				
				
			
		elif opt.kind == "out":
			for g in geneBED.bed_vals[chrom]:
				
				for t_index in xrange(len(testBED.bed_vals[chrom])):
					t = testBED.bed_vals[chrom][t_index];
					if conditions(t.start,t.end,g.start,g.end,overlap,fraction):

						break;

					elif t_index == len(testBED.bed_vals[chrom])-1:
			
						out =  chrom + "\t" + str(g.start)+"\t"+str(g.end) + "\t" + str(g.value) + "\n"; 
						total += 1;
						outputfile.write(out);

						
	outputfile.close();

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