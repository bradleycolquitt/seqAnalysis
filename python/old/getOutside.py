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

	line = infile.readline();
	line = line.strip();
	isoform_store = [line.split()];
	isoform_strand = 0;
	if isoform_store[0][5] == "+": isoform_strand=1;
	
	for line in infile:
	
		line = line.strip();
		sline = line.split();
		if sline[3] == isoform_store[0][3]:
			isoform_store.append(sline);
			
		else:
#			print "next"
#			print isoform_strand
			# If length of isoform store is greater than one,
			#	extract TSSs
			#	find most 5'
			#	write out isoform
			if len(isoform_store)>1:
#				print "multiple"
				test_list=[];
				ind = 0;
				
				if isoform_strand:
					for i in range(len(isoform_store)): test_list.append(isoform_store[i][1]);
					ind = test_list.index(min(test_list))
				else:
					for i in range(len(isoform_store)): test_list.append(isoform_store[i][2]);
					ind = test_list.index(max(test_list))
#					print "minus: ", isoform_store
#					sys.exit()
				out = "\t".join(isoform_store[ind]) + "\n";
							
			else:
#				print "single"
				out = "\t".join(isoform_store[0]) + "\n";
				
			outfile.write(out);
				
			isoform_store=[sline];
#			print "reset: ", isoform_store
			isoform_strand = 0;
			if isoform_store[0][5] == "+": isoform_strand=1;	
			
			


	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
