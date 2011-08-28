#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

"""
Take bed input with: chr start end strand value
Generate outfile with windows with user specified:
	Length
	Number
	Direction (upstream downstream)
	
	output format:
		chr start end strand value window-number
		
		--such that every input item will have N+1 output items
		--keep original interval as window-number 0
"""		

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	
	parser.add_option("-n", action="store", type="int", dest="number", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');

	for line in infile:
                
		line = line.strip();
		sline = line.split();
                new_end = atoi(sline[2]);
                new_start = atoi(sline[1]);
                length = int(round((new_end-new_start)/opt.number))
                #out = sline[0] + "\t" + sline[1] + "\t" + sline[2] + "\t" + sline[3] + "\t" + str(0) + "\t" + sline[5] + "\n"; 
                #outfile.write(out);
		for i in xrange(opt.number):
                        
                        out = [];
                        
                        if sline[5] == "+":
                                new_end = new_start+(length-1)
                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(i+1) + "\t" + sline[5] + "\n"
                                new_start = new_start+(length-1)
                        else:
                                new_start = new_end-(length-1)
                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(i+1) + "\t" + sline[5] + "\n"
                                new_end = new_end-(length-1);
                                             
                        outfile.write(out);
                                

	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
