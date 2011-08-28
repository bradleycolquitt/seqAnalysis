#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator





def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")

	(opt, args) = parser.parse_args(argv)

	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	
	par = 0
	antipar = 0
	total = 0
	
	# Check if strand represented as 1 or 0
	# If not, convert
	newin = []
	print "Formatting strands..."
	for line in infile:
		line = line.strip();
		sline = line.split();
		if type(sline[3]) != type(1):
			index_list = [3,6];
			for index in index_list:
				if sline[index] == "+":
					sline[index] = str(1);
				elif sline[index] == "-":
					sline[index] = str(0);
			newin.append("\t".join(sline)+"\n");
		infile = newin;
			
		
	
	print "Comparing strands..."
	for line in infile:
		line = line.strip();
		sline = line.split();
		if atoi(sline[3]) ^ atoi(sline[6]):
			antipar += 1
			total +=1
		else:
			par += 1
			total +=1
	
	output = "Total reads: " + str(total) + "\n" + "Parallel reads: " + str(par) + "\n" + "Antiparallel reads: " + str(antipar);
	outfile.write(output);

if __name__ == "__main__":
	main(sys.argv) 	
