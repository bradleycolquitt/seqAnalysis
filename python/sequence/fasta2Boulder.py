#!/usr/bin/env python

import sys, re, pdb
import Bio.SeqIO as seq

class boulder_io:
    def __init__(self, id, sequence):
        self.id = id
        self.template = sequence
        self.task = "generic"
        self.pick_left = 1
        self.pick_internal = 0
        self.pick_right = 1
        self.primer_opt_size = 20
        self.primer_min_size = 15
        self.primer_max_size = 25
        self.primer_max_n = 0
        seq_len = len(sequence)
        self.product_size = [seq_len - 80, seq_len - 20]
        
    def format(self):
        out = "\n".join(["SEQUENCE_ID={0}",
        "SEQUENCE_TEMPLATE={1}",
        "PRIMER_TASK={2}",
        "PRIMER_PICK_LEFT_PRIMER={3}",
        "PRIMER_PICK_INTERNAL_OLIGO={4}",
        "PRIMER_PICK_RIGHT_PRIMER={5}",
        "PRIMER_OPT_SIZE={6}",
        "PRIMER_MIN_SIZE={7}",
        "PRIMER_MAX_SIZE={8}",
        "PRIMER_MAX_NS_ACCEPTED={9}",
        "PRIMER_PRODUCT_SIZE_RANGE={10}-{11}",
        "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/usr/local/bin/primer3_config/",
        "=\n"])
        
        return(out.format(self.id,
                    self.template,
                    self.task,
                    self.pick_left,
                    self.pick_internal,
                    self.pick_right,
                    self.primer_opt_size,
                    self.primer_min_size,
                    self.primer_max_size,
                    self.primer_max_n,
                    self.product_size[0],
                    self.product_size[1]))
        

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1].split(".fa")[0] + ".boul", 'w')
    
    records = list(seq.parse(infile, "fasta"))
    new_seq = ""
    
    for record in records:
        boul = boulder_io(record.id, record.seq)
        outfile.write(boul.format()) 

if __name__ == '__main__':
    main(sys.argv)