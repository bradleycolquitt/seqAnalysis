#!/usr/bin/env python

import re, os, shutil, time, sys
from math import *
from string import *
from optparse import OptionParser
import operator


def main(argv):
	parser = OptionParser()
	
	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-t", "--head", action="store", type="string", dest="head", metavar="<str>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input);
	outputfile = open(opt.output, 'w');
	outputfile.write(">" + opt.head + "\n");
	for line in infile:
		if not re.match(">", line):
			line = line.strip();
			#sline = line.split();
			#print len(sline[0]);
			#print sline;
			#sys.exit(1);
			out = line;
			outputfile.write(out);

	outputfile.close();
		



if __name__ == "__main__":
	main(sys.argv)