	#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

"""
Three sets of intervals of length opt.window:
        opt.flank times upstream and downstream
        gene body divided into opt.number equally spaced windows of length opt.window
                
        
"""

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	sum = 0;

	for line in infile:
		line = line.strip();
		sline = line.split();
		temp = atoi(sline[2])-atoi(sline[1])
               	sum += temp;

	print sum;


if __name__ == "__main__":
	main(sys.argv) 	
