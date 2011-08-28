#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

def main(argv):

	parser = OptionParser()

	parser.add_option("-t", "--test", action="store", type="string", dest="test", metavar="<str>")
	parser.add_option("-g", "--gene", action="store", type="string", dest="gene", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<file>")


	(opt, args) = parser.parse_args(argv)

	test = open(opt.test, 'r');
	gene = open(opt.gene, 'r');
	outfile = open(opt.output, 'w');

	
	
	# Input two name2 files: one with coordinates, one without
	# Associate coordinates with file with none
	count = 0
	count_t = 0
	for line in test:
		count_t += 1;
		line = line.strip();
		sline = line.split();
		for gline in gene:
			gsline = gline.split();
			#print gsline[3]
			#print sline[0];
			#sys.exit(1);
			count += 1;
			if re.match(gsline[3], sline[0]):
				print "here"
				sys.exit(1);
				out = gsline[0] + "\t" + gsline[1] + "\t" + gsline[2] + "\t" + sline[1] + sline[2] + "\t" + sline[3] + "\n";
				outfile.write(out);
				#break;
	print count_t
	print count
	
if __name__ == "__main__":
	main(sys.argv) 	
