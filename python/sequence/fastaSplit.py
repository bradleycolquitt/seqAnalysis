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
	switch = 1;
	count = 0;
	for line in infile:
		line = line.strip();
		
		if switch:
			ssline = line.split("_");
			out = ssline[1] + "\t" + str(ssline[2]) + "\t" + str(ssline[3]) + "\t" + ssline[4] ;
			switch = 0;
		else:
			count += 1;
			line = line.strip();
			sline = line.split();
			out = str(count) + "\t" + sline[0] + "\t" + out + "\t" + sline[1] + "\n"; 
			print out;
			sys.exit(1);
			outputfile.write(out);
			switch = 1;

	outputfile.close();
		



if __name__ == "__main__":
	main(sys.argv)