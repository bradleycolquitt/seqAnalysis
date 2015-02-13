#! /usr/bin/env python

import sys, os

def main(argv):
    
    gtf_in = open(argv[1], 'r')
    gtf_name = os.path.basename(argv[1]).split(".gtf")[0]
    gtf_out = open(gtf_name + "_split", 'w')
    

    for line in gtf_in:
        line_tab = line.strip().split("\t")
        if len(line_tab) < 9:
            pass
        info = line_tab[8].split("; ")
        vals = []
        for field in info:
            vals.append(field.split()[1])
        
        out = "\t".join(line_tab[:8]) + "\t".join(vals) + "\n"
        gtf_out.write(out)

if __name__ == "__main__":
    main(sys.argv)
