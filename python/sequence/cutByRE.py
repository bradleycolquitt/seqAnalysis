#!/usr/bin/env python

import sys, re, pdb
import Bio.SeqIO as seq

def main(argv):
    
    infile = open(argv[1])
    reseq = argv[2]
    outfile = open(argv[1].split(".fa")[0] + "_split" + reseq + ".fa", 'a')
    
    records = list(seq.parse(infile, "fasta"))
    new_seq = ""
    
    for record in records:
        
        seq_split = record.seq.split(reseq)
        if len(seq_split) >= 3:
            for ind in range(1, len(seq_split)):
                if len(seq_split[ind]) >= 100:
                    new_seq = seq_split[ind]
        elif len(seq_split) == 2:
            #pdb.set_trace()
            continue
            new_seq = seq_split[0]
        elif len(seq_split) == 1:
            pdb.set_trace()
            continue
            new_seq = seq_split[0]
            
        record._set_seq(new_seq)
        #try:
        #    seq.write(record, outfile, "fasta")
        #except TypeError:
        #    pdb.set_trace()
    
    seq.write(records, outfile, "fasta")

if __name__ == '__main__':
    main(sys.argv)