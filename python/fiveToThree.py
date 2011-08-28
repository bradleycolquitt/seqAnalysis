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
			index_list = [3,6];
			for index in index_list:
				
				if sline[index] == "+":
					sline[index] = str(1);
				elif sline[index] == "-":
					sline[index] = str(0);
			newin.append("\t".join(sline)+"\n");
		infile = newin;
			
		

	
	pos = 0.0;
	orient = 0;
	for line in infile:
		
		line = line.strip();
		sline = line.split();

		readMid =(atoi(sline[5])+atoi(sline[4]))/2
		geneLength = atoi(sline[2])-atoi(sline[1]);
		
		#print sline;
		#print readMid
		#print geneLength
		#sys.exit(1);

		gene_pos = [atoi(sline[1])];
		start = atoi(sline[1]);
		
		for i in xrange(geneLength-1):
			start += 1;
			gene_pos.append(start);

		#print len(gene_pos);
		#sys.exit(1);
		for index in gene_pos:
			if index == readMid:
				#print index;
				#sys.exit(1);
				pos = (float(index)-atoi(sline[1]))/geneLength;
				#print pos;
				#print sline[1]
				#print sline[2]	
				#sys.exit(1);
				break;
		if atoi(sline[3]) ^ atoi(sline[6]):
			orient = 0;
		else:
			orient = 1; 
		
		out = sline[0] + "\t" + sline[1] + "\t" + sline[2] + "\t" + str(pos) + "\t" + str(orient) + "\n" 
 		
		outfile.write(out);

	outfile.close();			
			
			
		
		
if __name__ == "__main__":
	main(sys.argv) 	
