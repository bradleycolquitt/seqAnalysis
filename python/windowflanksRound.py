	#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from math import *
from optparse import OptionParser
import operator

"""
Three sets of intervals of length opt.window:
        opt.flank times upstream and downstream
        gene body divided into opt.number equally spaced windows of length opt.window
                
        
"""

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-w", action="store", type="int", dest="window", metavar="<int>")
	parser.add_option("-n", action="store", type="int", dest="number", metavar="<int>")
	parser.add_option("-d", action="store", type="string", dest="direction", metavar="<str>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	window = opt.window
	round_factor = int(-log10(window))
	number = opt.number
	direction = opt.direction
	flank_range = xrange(number*2)
	genome_sizes = open('/home/user/packages/bedtools/genomes/mouse.mm9.genome','r');

	d = {}
	for line in genome_sizes:
		line = line.strip();
		sline = line.split();
		if len(sline) > 0:
			d[sline[0]]=sline[1];
	
	for line in infile:
		line = line.strip();
		sline = line.split();
		strand = sline[5];
		initial_start = 0;
		initial_end = 0;
		
		if direction == "start":
			initial_start = int(round(atoi(sline[1]),round_factor)) + 1
			initial_end = int(round(atoi(sline[2]),round_factor))
			#initial_start = atoi(sline[1]) + 1
			#initial_end = atoi(sline[2])
		elif direction == "end":
			initial_end = int(round(atoi(sline[1]),round_factor)) 
			initial_start = int(round(atoi(sline[2]),round_factor)) + 1
			#initial_end = atoi(sline[1])
			#intial_start = atoi(sline[2]) + 1
		if strand == "+":
			start = initial_start - number * window
			for index in flank_range:
				end = start + window - 1
				out = sline[0] + "\t" + str(start) + "\t" + str(end) + "\t" + sline[3] + "\t" + str(index+1) + "\t" + sline[5] + "\n"
				if (start >= 0) and (end <= atoi(d[sline[0]])): 
                                        outfile.write(out)
				start = end + 1
				
		else:
			end = initial_end + number * window
			for index in flank_range:
				start = end - window + 1
				out = sline[0] + "\t" + str(start) + "\t" + str(end) + "\t" + sline[3] + "\t" + str(index+1) + "\t" + sline[5] + "\n"
				if (start >= 0) and (end <= atoi(d[sline[0]])): 
                                        outfile.write(out)
				end = start - 1
				

	outfile.close();
		
if __name__ == "__main__":
	main(sys.argv) 	
