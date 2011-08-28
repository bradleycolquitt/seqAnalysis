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

	
	


	print "Accessing test BED...";
	testBED = BED_mod.BED(species="mm9",file=opt.test,bed_type="BED_GRAPH",n=opt.name);
	print "Accessing gene BED...";
	geneBED = BED_mod.BED(species="mm9",file=opt.gene,bed_type="BED_GRAPH",n=opt.name);
	outputfile = open(opt.outFile, 'w');

	for chrom in testBED.keys():
		print "Test chromosome: " + chrom;
		if opt.kind == "in":
			for t in testBED.bed_vals[chrom]:
				for g in geneBED.bed_vals[chrom]:
					if conditions(t.start,t.end,g.start,g.end,overlap,fraction):
						out = chrom + "\t" + str(g.start)+"\t"+str(g.end)+"\t"+str(g.value) + "\t" + str(t.start)+ "\t" + str(t.end) + "\t" + str(t.value) + "\n";
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
						outputfile.write(out);

						
					
	print "Finishing..."
	outputfile.close();

def conditions(ts,te,gs,ge,o,f):
	a = 0;
	span = ge-gs;
	if ts<=gs and te>=ge:
		a = 1;
	elif ts>=gs and te<=ge:
		if (te-ts)/span>=f: 
			a=1;
	elif ts<=gs and te>=gs+o:
		if (te-gs)/span>=f:
			a=1;
	elif ts<=ge+o and te>=ge:
		if (ge-ts)/span>=f: 
			a=1;
			x=te-ts
			print ts,te,gs,ge,o,f,x,span, a
			sys.exit(1)
	
	
	
		
	return a;

			

if __name__ == "__main__":
	main(sys.argv)