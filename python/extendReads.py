#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-e", "--extend", action="store", type="int", dest="extend", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	strand = 3;

	sline = infile.readline().strip().split()
	for index in xrange(len(sline)):
		if sline[index]=="+" or sline[index]=="-":
			strand=index
			break;
		
	infile.close()
	infile = open(opt.input, 'r')
	for line in infile:
		line = line.strip();
		sline = line.split("\t");
		if sline[strand] == "+":
			sline[2] = str(atoi(sline[1])+opt.extend);
		elif sline[strand] == "-":
			sline[1] = str(atoi(sline[2])-opt.extend);	
		out = "\t".join(sline) + "\n";
		outfile.write(out);

	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
