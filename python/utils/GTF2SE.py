#!/usr/bin/env python

# Input: GTF file
# Output: GFF3 file
#   For each internal exon (example from MISO site)
#   chr  SE      gene    triple_left triple_right .       strand       .       ID=gene_id:exon_number;Name=gene_id:exon_number
#   chr  SE      mRNA    triple_left triple_right .       strand       .       ID=gene_id:exon_number.A;Parent=gene_id:exon_number
#   chr  SE      mRNA    triple_left triple_right .       strand       .       ID=gene_id:exon_number.B;Parent=gene_id:exon_number
#   chr  SE      exon    up_left up_right .       strand       .       ID=gene_id:exon_number.A.up;Parent=gene_id:exon_number.A
#   chr  SE      exon    mid_left mid_right .       strand       .       ID=gene_id:exon_number.A.se;Parent=gene_id:exon_number.A
#   chr  SE      exon    down_left down_right .       strand       .       ID=gene_id:exon_number.A.dn;Parent=gene_id:exon_number.A
#   chr  SE      exon    up_left up_right .       strand       .       ID=gene_id:exon_number.B.up;Parent=gene_id:exon_number.B
#   chr  SE      exon    down_left down_right .       strand       .       ID=gene_id:exon_number.B.dn;Parent=gene_id:exon_number.B

import sys
import pdb

def roll_list_back(l, new_val):
    for i in range(1,len(l)):
        l[i-1] = l[i]
    l[len(l) - 1] = new_val
    return(l)
    
class gtf_field:
    def __init__(self, chr, source, feature, start, end, strand, gene_id, transcript_id, exon_number):
        self.chr = chr
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.exon_number = exon_number
        
def parse_gtf(line):
    sline = line.split()
    return(gtf_field(sline[0], sline[1], sline[2], sline[3], sline[4], sline[6],
                     format_annotation(sline[9]), format_annotation(sline[11]),
                     format_annotation(sline[13])))

def format_annotation(anno):
    anno = anno.split(";")[0]
    #anno = anno.split("\"")[1]
    return(anno)
    
def out_format(exon_set):
    #pdb.set_trace()
    strand = exon_set[0].strand
    id = ":".join([exon_set[1].transcript_id, exon_set[1].exon_number])
    out = ["","","","","","","",""]
    exon_ind = ""
    if strand == "+":
        exon_ind = (0,1,2)
    elif strand == "-":
        exon_ind = (2,1,0)
        
    out[0] = "{0}\tSE\t{1}\t{2}\t{3}\t.\t{4}\t.\tID={5};Name={6}\n".format(exon_set[0].chr,
                                                                             "gene",
                                                                             exon_set[exon_ind[0]].start,
                                                                             exon_set[exon_ind[2]].end,
                                                                             strand,
                                                                             id,
                                                                             id)
    out[1] = "{0}\tSE\t{1}\t{2}\t{3}\t.\t{4}\t.\tID={5};Parent={6}\n".format(exon_set[0].chr,
                                                                             "mRNA",
                                                                             exon_set[exon_ind[0]].start,
                                                                             exon_set[exon_ind[2]].end,
                                                                             strand,
                                                                             id + ".A",
                                                                             id)
    out[2] = "{0}\tSE\t{1}\t{2}\t{3}\t.\t{4}\t.\tID={5};Parent={6}\n".format(exon_set[0].chr,
                                                                             "mRNA",
                                                                             exon_set[exon_ind[0]].start,
                                                                             exon_set[exon_ind[2]].end,
                                                                             strand,
                                                                             id + ".B",
                                                                             id)
    out[3] = "{0}\tSE\t{1}\t{2}\t{3}\t.\t{4}\t.\tID={5};Parent={6}\n".format(exon_set[0].chr,
                                                                             "exon",
                                                                             exon_set[exon_ind[0]].start,
                                                                             exon_set[exon_ind[0]].end,
                                                                             strand,
                                                                             id + ".A.up",
                                                                             id + ".A")
    out[4] = "{0}\tSE\t{1}\t{2}\t{3}\t.\t{4}\t.\tID={5};Parent={6}\n".format(exon_set[1].chr,
                                                                             "exon",
                                                                             exon_set[exon_ind[1]].start,
                                                                             exon_set[exon_ind[1]].end,
                                                                             strand,
                                                                             id + ".A.se",
                                                                             id + ".A")
    out[5] = "{0}\tSE\t{1}\t{2}\t{3}\t.\t{4}\t.\tID={5};Parent={6}\n".format(exon_set[2].chr,
                                                                             "exon",
                                                                             exon_set[exon_ind[2]].start,
                                                                             exon_set[exon_ind[2]].end,
                                                                             strand,
                                                                             id + ".A.dn",
                                                                             id + ".A")
    out[6] = "{0}\tSE\t{1}\t{2}\t{3}\t.\t{4}\t.\tID={5};Parent={6}\n".format(exon_set[0].chr,
                                                                             "exon",
                                                                             exon_set[exon_ind[0]].start,
                                                                             exon_set[exon_ind[0]].end,
                                                                             strand,
                                                                             id + ".B.up",
                                                                             id + ".B")
    out[7] = "{0}\tSE\t{1}\t{2}\t{3}\t.\t{4}\t.\tID={5};Parent={6}\n".format(exon_set[0].chr,
                                                                             "exon",
                                                                             exon_set[exon_ind[2]].start,
                                                                             exon_set[exon_ind[2]].end,
                                                                             strand,
                                                                             id + ".B.dn",
                                                                             id + ".B")
    return(out)
    
#http://stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical    
def checkEqual3(lst):
       return lst[1:] == lst[:-1]
       
def main(argv):
    infile = open(argv[1], 'r')
    outfile = open(argv[2], 'w')
    exon_set = [0,0,0]
    first = True
    for line in infile:
        curr_gtf = parse_gtf(line)
        exon_set = roll_list_back(exon_set, curr_gtf)
        
        if not exon_set[0] == 0:
            if not checkEqual3([exon_set[0].transcript_id, exon_set[1].transcript_id, exon_set[2].transcript_id]):
                exon_set[0] = 0
            else:
                #if exon_set[1].transcript_id == "ENSMUST00000017981":
                #    pdb.set_trace()
                formatted_list = out_format(exon_set)
                for field in formatted_list:
                    outfile.write(field)
                    
        
if __name__ == "__main__":
    main(sys.argv)