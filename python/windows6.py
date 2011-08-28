	#!/usr/bin/env python

import re, os, shutil, time, sys
from string import *
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
	#parser.add_option("-f", action="store", type="int", dest="flank", metavar="<int>")
	parser.add_option("-w", action="store", type="int", dest="window", metavar="<int>")
	#parser.add_option("-n", action="store", type="int", dest="number", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	
	infile = open(opt.input, 'r');
	outfile = open(opt.output, 'w');
	window = opt.window
	#number = opt.number
	#flank = opt.flank
	#fsect = xrange(flank)
	#gsect = xrange(number)
	

	for line in infile:
               
		line = line.strip();
		sline = line.split();
                new_end = atoi(sline[2]);
                new_start = atoi(sline[1]);
		advance = xrange((new_end-new_start+1)/window)
                #remain = (atoi(sline[2])-atoi(sline[1])) - number*window
                #if remain < 0:
                #        continue;
                #rwindow = remain/(number-1)
                #print sline, remain, rwindow
                #sys.exit()
        
##                if opt.direct == "up":
##                        start = atoi(sline[1])
##                        end = atoi(sline[2])
##                elif opt.direct == "down":
##                        start = atoi(sline[2])
##                        end = atoi(sline[1])
                       
                """
                for x in xrange(3):
                        if x == 0:
                                new_start = new_start-flank*window
                                new_end = new_end+flank*window
                                for i in fsect:
                                        out = "";
                                        if sline[5] == "+":
                                                new_end = new_start+(window-1)
                                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(i+1) + "\t" + sline[5] + "\n"
                                                new_start = new_end+1
                                                
                                                
                                        else:
                                                new_start = new_end-window+1
                                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(i+1) + "\t" + sline[5] + "\n"
                                                new_end = new_start-1;   
                                        outfile.write(out);
                        elif x == 1:
		"""
               	new_start = atoi(sline[1])
                new_end = atoi(sline[2])
                for n in advance:
               		out = "";
                        if sline[5] == "+":
                        	new_end = new_start + window-1
                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(n+1) + "\t" + sline[5] + "\n"
                                new_start = new_end + 1
                        else:
                                new_start = new_end - window+1
                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(picn+1) + "\t" + sline[5] + "\n"
                                new_end = new_start - 1
                        outfile.write(out);
                """                               
                        elif x == 2:
                                new_start = atoi(sline[2])
                                new_end = atoi(sline[1])
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
                                        outfile.write(out);
##                                 for i in sect:        
##                                        out = "";
##                                        if sline[5] == "+":
##                                                new_end = new_start+(window-1)
##                                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(abs(i-flank)) + "u" + "\t" + sline[5] + "\n"
##                                                new_start = new_end+1
##                                                
##                                                
##                                        else:
##                                                new_start = new_end-window+1
##                                                out = sline[0] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + sline[3] + "\t" + str(abs(i-flank)) + "u" + "\t" + sline[5] + "\n"
##                                                new_end = new_start-1;   
##                                        outfile.write(out);
		 """

	outfile.close();

if __name__ == "__main__":
	main(sys.argv) 	
