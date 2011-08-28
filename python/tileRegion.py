#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-w", action="store", type="int", dest="window", metavar="<int>")
	parser.add_option("-s", action="store", type="int", dest="step", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	window = opt.window
	step = opt.step
	
	for line in infile:
                
		line = line.strip();
		sline = line.split();
		new_end = atoi(sline[1])+window-1;
		new_start = atoi(sline[1]);
               
		while new_end <= atoi(sline[2]):
					
					out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\n";
					outfile.write(out);
					new_start = new_start+step
					new_end = new_start+(window-1)
					
                                

	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
