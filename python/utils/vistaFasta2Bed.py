#!/usr/bin/env python

import sys, re, pdb

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1].split(".fa")[0] + ".bed", 'w')
    
    expr = ""
   # expr2 = []
    species = ""
    element = ""
    for line in infile:
        if re.search(">", line):
            
            line = line.strip().split("|")
            #pdb.set_trace()
            if line[0] == ">Human":
                species = "hs"
            elif line[0] == ">Mouse":
                species = "mm"
                
            
            element = line[2].split(" ")
            element = species + element[2]
                
            #if re.search("positive", line[3]): 
            #    expr = line[4:]
            #    for tissue in expr:
            #        tissue = tissue.split("[")
            #        expr2.append(tissue[0].strip())
            #else:
            #    expr2 = ["negative"]
            
            pos = line[1].split(":")
            chr = pos[0]
            pos2 = pos[1].split("-")
            start = pos2[0]
            end = pos2[1]
            
            outline = "\t".join([chr, start, end, element, "0", "+"]) + "\n"
            outfile.write(outline)
                
  #          expr2 = []

if __name__ == '__main__':
    main(sys.argv)