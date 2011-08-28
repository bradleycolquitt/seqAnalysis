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
	for line in infile:
		line = line.strip();
		sline = line.split();
		if type(sline[3]) != type(1):
			index_list = [3,4];
			for index in index_list:
				
				if sline[index] == "+":
					sline[index] = str(1);
				elif sline[index] == "-":
					sline[index] = str(0);
			newin.append("\t".join(sline)+"\n");
		infile = newin;
			
		

	
	par = 0.0;
	subtotal = 0.0;
	prev = "1st";
	oline=[1,2];
	score = 0.0;
	for line in infile:
		
		line = line.strip();
		sline = line.split();
		
		if oline[1] == sline[1]:
			if atoi(sline[3]) ^ atoi(sline[4]):
				subtotal += 1
			else:
				par += 1
				subtotal += 1
			oline = sline;
		elif oline[1] == 2:
			if atoi(sline[3]) ^ atoi(sline[4]):
				subtotal += 1
			else:
				par += 1
				subtotal += 1
			oline = sline;

		else:
			if subtotal > 0:
				score = float(par)/float(subtotal);
				output = str(oline[0]) + "\t" + str(oline[1]) + "\t" + str(oline[2]) + "\t" + str(score) + "\t" + str(subtotal) + "\n"; 	
				outfile.write(output);
			
			par = 0.0;
			subtotal = 0.0;
			if atoi(sline[3]) ^ atoi(sline[4]):
				subtotal += 1
			else:
				par += 1
				subtotal += 1
			oline = sline;
	if subtotal > 0:
		score = float(par)/float(subtotal);
		output = str(oline[0]) + "\t" + str(oline[1]) + "\t" + str(oline[2]) + "\t" + str(score) + "\t" + str(subtotal) + "\n"; 	
		outfile.write(output);

			
			
			
		
		
if __name__ == "__main__":
	main(sys.argv) 	
