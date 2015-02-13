#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	#parser.add_option("-d", action="store", type="str", dest="direct", metavar="<str>")
	parser.add_option("-l", action="store", type="int", dest="length", metavar="<int>")
	#parser.add_option("-n", action="store", type="int", dest="number", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	length = opt.length
	qual = "I"*length;		
	
	for line in infile:
		line = line.strip();
		sline = line.split();
		out = "@test" + "\n" + sline[1] + "\n" + "+" + "\n" + qual + "\n"
		outfile.write(out);
			
	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
