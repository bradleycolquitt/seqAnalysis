	#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

"""
Format set of BED regions for index at given resolution.
        Read in BED file with chr, start, stop, name, score, strand
        Round start and stop to nearest window resolution
        Extend from start or stop (depending on strand and specified position) a given number of windows
        
"""

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-d", action="store", type="int", dest="flank", metavar="<int>")
	parser.add_option("-w", action="store", type="int", dest="window", metavar="<int>")
	parser.add_option("-n", action="store", type="int", dest="number", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	
	window = opt.window
	round_factor = -log(window,10)
        number = opt.number
	nsect = xrange(number)
	genome_sizes = open('/users/bradcolquitt/BEDTools/genomes/mouse.mm9.genome.txt','r');

	d = {}
	for line in genome_sizes:
		line = line.strip();
		sline = line.split();
		if len(sline) > 0:
			d[sline[0]]=sline[1];
	

	for line in infile:
		line = line.strip();
		sline = line.split();
		
		new_start = round(atoi(sline[1]),-2) + 1;
                new_end = round(atoi(sline[2]),-2) + 1;
		out = "";
                for n in gsect:	
			if sline[5] == "+":
				new_end = new_start + window-1
				out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(flank+n+1) + "\t" + sline[5] + "\n"
				new_start = new_end + rwindow + 1
			else:
				new_start = new_end - window+1
				out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(flank+n+1) + "\t" + sline[5] + "\n"
				new_end = new_start - rwindow - 1
			if (new_start and new_end >=0) and (new_start and new_end <= atoi(d[sline[0]])): 
				outfile.write(out);
                

	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
