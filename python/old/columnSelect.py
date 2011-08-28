#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator





def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-c", "--columns", action="store", type="string", dest="columns", metavar="<str>")
#	parser.add_option("-fc", "--filler_columns", action="store", type="string", dest="filler_columns", metavar="<str>")
#	parser.add_option("-f", "--filler", action="store", type="string", dest="filler", metavar="<str>")

	(opt, args) = parser.parse_args(argv)

	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	col = opt.columns.split(',');
#	filler_columns = opt.filler_columns.split(',')
#	filler = opt.columns.split(',');
	
	for line in infile:
		line = line.strip();
		sline = line.split();
		items = [];
		for c in col:
			index = atoi(c)-1;
			items.append(sline[index]);
		out = "\t".join(items) + "\n";
		outfile.write(out);
		

	outfile.close();

if __name__ == "__main__":
	main(sys.argv)
		