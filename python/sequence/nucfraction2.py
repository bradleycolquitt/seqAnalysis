#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator





def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-p", "--pattern", action="store", type="string", dest="pattern", metavar="<str>")


	(opt, args) = parser.parse_args(argv)

	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	pl = list(opt.pattern.lower())
	
	pu = list(opt.pattern.upper())
	
	pa = pl + pu
	
	
	# Read in fasta, split, take sequence (element 1)
	# split pattern
	# iterate across sequence
	# if pattern elements or, +1
	# output 

	for line in infile:
		count = 0;
		line = line.strip();
		sline = line.split();
		for p in pa:
			count += sline[6].count(p);
		
		fraction = float(count)/len(sline[6])	
		
		out = "\t".join(sline[0:6])  + "\t" + str(fraction) + "\n";
		outfile.write(out);
		
	
if __name__ == "__main__":
	main(sys.argv) 	
