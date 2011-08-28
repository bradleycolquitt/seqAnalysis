#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

 

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o",action="store",type="string",dest="output",metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	
	chrom = "";
	ind = 1;
	span = "";
	for line in infile:
			
		line = line.strip();
		sline = line.split();
		if ind:
			span = str(atoi(sline[2])-atoi(sline[1])+1)
			ind = 0;
		if sline[0] == chrom:
			out = str(atoi(sline[1])+1) + "\t" + sline[3] + "\n"
			outfile.write(out);
		else:	
			out = "variableStep" + "\t" + "chrom=" + sline[0] + "\t" + "span=" + span + "\n" + str(atoi(sline[1])+1) + "\t" + sline[3] + "\n"
			outfile.write(out);
			chrom = sline[0];
	
	outfile.close();
	


	


if __name__ == "__main__":
	main(sys.argv) 	
