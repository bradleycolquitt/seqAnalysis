#!/usr/bin/env python

import re, os, shutil, time, sys
from math import *
from string import *
from optparse import OptionParser
import operator


def main(argv):
	parser = OptionParser()
	
	parser.add_option("-a", action="store", type="string", dest="input1", metavar="<str>")
	parser.add_option("-b", action="store", type="string", dest="input2", metavar="<str>")
	parser.add_option("-o", action="store", type="string", dest="output", metavar="<str>")

	(opt, args) = parser.parse_args(argv)
	
	infile1 = open(opt.input1, 'r');
	infile2 = open(opt.input2, 'r');
	outputfile = open(opt.output, 'w');
	for i in xrange(len(infile1)):
		out = infile1.readline() + "\t" + infile2.readline();
		outputfile.write(out);
	outputfile.close();
		



if __name__ == "__main__":
	main(sys.argv)