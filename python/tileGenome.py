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
	parser.add_option("-w", action="store", type="int", dest="window", metavar="<int>")
	#parser.add_option("-n", action="store", type="int", dest="number", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	window = opt.window
	#number = opt.number
	#sect = xrange(number*2)
	
	d = {}
	for line in infile:
		line = line.strip();
		sline = line.split();
		if len(sline) > 0:
			d[sline[0]]=sline[1];
		
	
	for c in d.keys():
		start = 0;
		end = window-1;
		ad = window/2;
		term = atoi(d[c]);
		switch = 1;
		while switch:
			if end > term: switch = 0;
			else:
				out = c + "\t" + str(start) + "\t" + str(end) + "\n"
				outfile.write(out);
				start = start + ad;
				end = start + window - 1;
	                                

	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
