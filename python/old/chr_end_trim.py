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

	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<str>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>")

	
	(opt, args) = parser.parse_args(argv)
	type = "BED_GRAPH";
	infile = open(opt.infile, 'r');
	for line in infile:
		line = line.strip();
		sline = line.split();
		if len(sline) == 3:
			type = "BED3";
			break;
		else:
			break;
	print "Accessing BED...";
	inBED = BED_mod.BED(species="mm9",file=opt.infile,bed_type=type);
	outputfile = open(opt.outfile, 'w');	
	
	if type == "BED_GRAPH":
		out = "track type=\"bedGraph\" name=\"" + str(opt.infile) + "\" " + "description=\"" + str(opt.outfile) + "\"\n"  ;
		outputfile.write(out);
	elif type =="BED3":
		out = "track name=\"" + str(opt.infile) + "\"\t" + "description=\"" + str(opt.outfile) + "\"\n"  ;	
		outputfile.write(out);
	for chrom in inBED.keys():
		for t in xrange(len(inBED.bed_vals[chrom])-1):
			if type == "BED_GRAPH":
				out = chrom + "\t" + str(inBED.bed_vals[chrom][t].start) + "\t" + str(inBED.bed_vals[chrom][t].end) + "\t" + str(inBED.bed_vals[chrom][t].value) + "\n"  ;
				outputfile.write(out);

			elif type == "BED3":
				out = chrom + "\t" + str(inBED.bed_vals[chrom][t].start) + "\t" + str(inBED.bed_vals[chrom][t].end) + "\n" ;
				outputfile.write(out);
				

	outputfile.close()
				

if __name__ == "__main__":
	main(sys.argv)