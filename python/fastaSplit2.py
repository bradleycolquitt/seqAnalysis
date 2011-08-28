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

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outputfile = open(opt.output, 'w');

	for line in infile:
		line = line.strip();
		sline = line.split();
		ssline = sline[0].split("_");
		 
		out = ssline[1] + "\t" + ssline[2] + "\t" + ssline[3] + "\t" + ssline[4] + "\t" + "\t".join(sline[1:]) + "\n";
		outputfile.write(out);
	outputfile.close();
		



if __name__ == "__main__":
	main(sys.argv)
