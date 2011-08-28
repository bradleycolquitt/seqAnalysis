#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator





def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	#parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-c", "--columns", action="store", type="string", dest="columns", metavar="<str>")

	(opt, args) = parser.parse_args(argv)

	infile = open(opt.input, 'r');
	#outfile = open(opt.output, 'w');
	col = opt.columns
	items=[]
	for line in infile:
		line = line.strip();
		sline = line.split();
		index = atoi(col)-1;
		items.append(atof(sline[index]));
	avg = sum(items)/len(items)
	print sum(items)
	print len(items)
	print avg;			

	

if __name__ == "__main__":
	main(sys.argv)
		