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
	
	
	# Read in fasta, split, take sequence (element 1)
	# iterate across sequence combining index, index+1
	# check if combination == "CG", +1 to count if true
	# output 

	for line in infile:
		line = line.strip();
		sline = line.split();
		ssline = sline[1].split();
		print len(ssline);
		sys.exit(1);
		if type(sline[3]) != type(1):
			index_list = [3,4];
			for index in index_list:
				if sline[index] == "+":
					sline[index] = str(1);
				elif sline[index] == "-":
					sline[index] = str(0);
			newin.append("\t".join(sline)+"\n");
		infile = newin;
			
		
	

	for line in infile:
		antipar = 0;
		par = 0;
		subtotal = 0;

		line = line.strip();
		sline = line.split();

		if atoi(sline[3]) ^ atoi(sline[4]):
			antipar += 1
			subtotal += 1
		else:
			par += 1
			subtotal += 1
		for line2 in infile:
			line2 = line2.strip();
			sline2 = line2.split();
			if sline[1] == sline2[1]:
				if atoi(sline2[3]) ^ atoi(sline2[4]):
					antipar += 1
					subtotal +=1
				else:
					par += 1
					subtotal +=1
			else:
				break;
		score = par/subtotal
		output = str(line[0]) + "\t" + str(line[1]) + "\t" + str(score) + "\n"; 	
		outfile.write(output);

if __name__ == "__main__":
	main(sys.argv) 	
