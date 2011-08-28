#! /usr/bin/env python

"""
For each item in query file, find closest item in subject
Direction indicates whether to look upstream, downstream, or both
Input files should be bed4 or bed6
"""

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', dest="query")
    parser.add_argument('-s', dest="subject")
    parser.add_argument('-d', dest="direct")
    args = parser.parse_args()
    query = open(args.query)
    subject = open(args.subject)
    
    ## Split subject by chrom, build dict
    subject_dict = {}
    for line in subject:
        sline = line.strip().split()
        if sline[0] not in subject_dict.keys():
            subject_dict[sline[0]] = [sline]
        else:
            subject_dict[sline[0]] = subject_dict[sline[0]].append(sline)
    
    for line in query:
        line = line.strip().split()
        
    

if __name__ == "__main__":
    main(sys.argv)
