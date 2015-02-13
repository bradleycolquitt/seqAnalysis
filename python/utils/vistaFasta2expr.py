#!/usr/bin/env python

import sys, re, pdb

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1].split(".fa")[0] + "_expr.txt", 'w')
    
    expr = ""
    expr2 = []
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
                
            if re.search("positive", line[3]): 
                expr = line[4:]
                for tissue in expr:
                    tissue = tissue.split("[")
                    expr2.append(tissue[0].strip())
            else:
                expr2 = ["negative"]
                
            outline = "\t".join([element, "-".join(expr2)]) + "\n"
            outfile.write(outline)
                
            expr2 = []

if __name__ == '__main__':
    main(sys.argv)