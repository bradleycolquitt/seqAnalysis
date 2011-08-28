#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
from optparse import OptionParser
import operator

def main(argv):

	parser = OptionParser()

	parser.add_option("-i", "--input", action="store", type="string", dest="input", metavar="<str>")
	parser.add_option("-o", "--ouput", action="store", type="string", dest="output", metavar="<str>")
	parser.add_option("-d", action="store", type="str", dest="direct", metavar="<str>")
	parser.add_option("-w", action="store", type="int", dest="window", metavar="<int>")
	parser.add_option("-n", action="store", type="int", dest="number", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	window = opt.window
	number = opt.number
	sect = xrange(number*2)
	

	for line in infile:
                
		line = line.strip();
		sline = line.split();
                new_end = atoi(sline[1]);
                new_start = atoi(sline[2]);
                
        
                if opt.direct == "up":
                        start = atoi(sline[1])
                        end = atoi(sline[2])
                elif opt.direct == "down":
                        start = atoi(sline[2])
                        end = atoi(sline[1])
                       
                new_start = start-number*window
                new_end = end+number*window-1        
                for i in sect:
                        out = [];
                
                        if sline[5] == "+":
                                new_end = new_start+(window-1)
                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(i-number+1) + "\t" + sline[5] + "\n"
                                new_start = new_end+1
                                
                        else:
                                new_start = new_end-window+1
                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(i-number+1) + "\t" + sline[5] + "\n"
                                new_end = new_start-1;

                                   
                        outfile.write(out);
                                

	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
