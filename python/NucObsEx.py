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
	pattern = re.compile(opt.pattern,re.I);
	
	# Read in fasta, split, take sequence (element 1)
	# split pattern
	# iterate across sequence
	# if pattern elements or, +1
	# output 

	for line in infile:
		fraction = 0;
		count1 = 0;
		count2 = 0;
		line = line.strip();
		sline = line.split();
	
		count1 += sline[4].count(pa[0]);
		count1 += sline[4].count(pa[2]);
		count2 += sline[4].count(pa[1]);
		count2 += sline[4].count(pa[3]);
		
		exp = count1 * count2
	
		m = pattern.findall(sline[4])
		obs = len(m)
		
		if exp > 0:
			fraction = (float(obs)/exp)*len(sline[4])	
		
		out = "\t".join(sline[0:4]) + "\t" + "\t".join(sline[5:]) + "\t" + str(fraction) + "\n";
		outfile.write(out);
		
	
if __name__ == "__main__":
	main(sys.argv) 	
