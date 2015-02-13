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
	parser.add_option("-f", action="store", type="int", dest="flank", metavar="<int>")
	parser.add_option("-w", action="store", type="int", dest="window", metavar="<int>")
	parser.add_option("-n", action="store", type="int", dest="number", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	window = opt.window
	round_factor = int(-log10(window))
	number = opt.number
	flank = opt.flank
	fsect = xrange(flank)
	gsect = xrange(number)
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
                initial_start = int(round(atoi(sline[1]),round_factor)) + 1;
		initial_end = int(round(atoi(sline[2]),round_factor));
                remain = (atoi(sline[2])-atoi(sline[1])) - number*window
                if remain < 0:
                        continue;
                rwindow = remain/(number-1)
        
                for x in xrange(3):
                        if x == 0:
                                flank_start_ = initial_start-flank*window
                                flank_end = initial_end+flank*window
                                for i in fsect:
                                        out = "";
                                        if sline[5] == "+":
                                                flank_end = flank_start+(window-1)
                                                out = sline[0] + "\t" + str(flank_start) + "\t" + str(flank_end) + "\t" + sline[3] + "\t" + str(i+1) + "\t" + sline[5] + "\n"
                                                flank_start = flank_end+1
                                                
                                                
                                        else:
                                                flank_start = flank_end-window+1
                                                out = sline[0] + "\t" + str(flank_start) + "\t" + str(flank_end) + "\t" + sline[3] + "\t" + str(i+1) + "\t" + sline[5] + "\n"
                                                flank_end = flank_start-1;  
					if (flank_start and flank_end >=0) and (flank_start and flank_end <= atoi(d[sline[0]])): 
                                        	outfile.write(out);
                        elif x == 1:
				new_start = initial_start
				new_end = initial_end
                                for n in gsect:
                                        out = "";
                                        if sline[5] == "+":
                                                new_end = new_start + window-1
                                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(flank+n+1) + "\t" + sline[5] + "\n"
                                                new_start = int(round(new_end + rwindow + 1,round_factor))
                                        else:
                                                new_start = new_end - window+1
                                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(flank+n+1) + "\t" + sline[5] + "\n"
                                                new_end = int(round(new_start - rwindow - 1))
                                        if (new_start and new_end >=0) and (new_start and new_end <= atoi(d[sline[0]])): 
                                        	outfile.write(out);
                                                
                        elif x == 2:
				
                                for i in fsect:
                                        out = "";
                                        if sline[5] == "+":
                                                new_end = new_start+(window-1)
                                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(flank+number+i+1)  + "\t" + sline[5] + "\n"
                                                new_start = new_end+1
                                        else:
                                                new_start = new_end-window+1
                                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(flank+number+i+1)  + "\t" + sline[5] + "\n"
                                                new_end = new_start-1;   
                                        if (new_start and new_end >=0) and (new_start and new_end <= atoi(d[sline[0]])): 
                                        	outfile.write(out);

	outfile.close();
		
if __name__ == "__main__":
	main(sys.argv) 	
