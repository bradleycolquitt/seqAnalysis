#! /usr/bin/env python

import os, sys
import pysam
import bam2bed
import sam
import argparse
import BAMtoWIG
import pdb
from subprocess import Popen

fastq_dir = "/media/storage3/data/fastq"
bam_dir = "/media/storage2/data/rna/tophat"
#bam_dir = "/home/user/data/rna/tophat"
wig_dir = "/media/storage2/data/wig"

class tophat_class:
    def __init__(self, date, sample, single_end, mean, sd, gtf, library_type, species):
        self.date = date
        self.sample = sample
        self.single_end = single_end
        self.input1 = ""
        self.input2 = ""
        if single_end:
            self.input1 = "/".join([fastq_dir, date, sample, index[0],
                                    "_".join([self.sample, '1.fastq'])])
        else:
            self.input1 = "/".join([fastq_dir, date, sample, 
                                  "_".join([self.sample, 'r1.fastq.gz'])])
            self.input2 = "/".join([fastq_dir, date, sample, 
                                    "_".join([self.sample, 'r2.fastq.gz'])])
        self.mean = mean
        self.sd = sd
        self.gtf = gtf
        self.library_type = library_type
        self.species = species
        self.to_map = True
        self.output = "/".join([bam_dir, self.date, self.sample])
        if os.path.exists(self.output):
            dec = raw_input("Output file exists. Bypass mapping [y/n]? ")
            if dec == 'y': self.to_map = False
        else: os.makedirs(self.output)
        self.log = "/".join([self.output, "log"])
        self.errorlog = "/".join([self.output, "error_log"])
        
    def map(self):        
        if self.to_map:
            if self.single_end:
                if self.gtf != 'blank': 
                    cmd_args = ['tophat', '-p', '6', '-o', self.output, 
                        '-r', self.mean, '--mate-std-dev', self.sd, 
                        '-G', self.gtf, '--library-type', self.library_type, 
                        self.species, self.input1]
                else:
                    cmd_args = ['tophat', '-p', '6', '-o', self.output, 
                        '-r', self.mean, '--mate-std-dev', self.sd,
                        '--library-type', self.library_type,
                        self.species, self.input1]
            else:
                if self.gtf != 'blank': 
                    cmd_args = ['tophat2',
                                '-p', '6',
                                '-o', self.output, 
                                '-r', self.mean,
                                '--mate-std-dev', self.sd,
                                '--max-multihits', '1',
                                '-G', '/home/user/lib/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/genes_chr.gtf',
                                '--transcriptome-index', '/home/user/lib/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/',
                                #'--prefilter-multihits',
                                '--no-coverage-search',
                                '--library-type', self.library_type, 
                                self.species, self.input1, self.input2]
                else:
                    cmd_args = ['tophat2',
                                '-o', self.output,
                                '-p', '6',
                                '-r', self.mean,
                                '--mate-std-dev', self.sd,
                                '--library-type', self.library_type,
                                self.species, self.input1, self.input2]
            
            log = open(self.log, 'a')
            #errorlog = open(self.errorlog, 'a')
            print "Mapping with tophat: " + " ".join(cmd_args[1:])
            print >>log, "Mapping with tophat: " + " ".join(cmd_args[1:])
            tophat = Popen(cmd_args, stdout=log, stderr=log)
            tophat.wait()
            #errorlog.close()
            #sam.proc(self.output + "/accepted_hits.bam", sort=False, rmdup="False", errorlog=errorlog)
            
    def wig(self):
        #pdb.set_trace()
        window_size = str(25)
        bamname = "/".join([self.output, "accepted_hits.bam"])
        wigdir = "/".join([wig_dir, self.sample])
        if not os.path.exists(wigdir): os.makedirs(wigdir)
        wigfile = "/".join([wigdir, "_".join([self.sample, window_size]) + ".wig"])
        extend = self.mean
        pe = not self.single_end
        wi = BAMtoWIG.windower(bamname, wigfile, window_size, extend, pe, False, False)
        print "Writing WIG..."
        wi.window()
        wi.wigfile.close()
        wi.tdf()
        
def tophat(date, sample, single_end, mean, sd, gtf, library_type, species):
    tophat_obj = tophat_class(date, sample, single_end, mean, sd, gtf, library_type, species)
    tophat_obj.map()
    
def main(argv):
     
    parser = argparse.ArgumentParser(description="Run tophat on given fastq files.\nOutputs to " + bam_dir)
    parser.add_argument('-d', dest='date')
    parser.add_argument('-i', '--input_fastq', nargs='+', metavar='fastq', required=True, 
                        dest='fastq_files', help='read1 and read3 of library')
    
    parser.add_argument('-r', '--mate-inner-dist', metavar='mean', required=True, 
                        dest='mean', help='mean size of library (minus adaptors)')
    parser.add_argument('-s', '--mate-std-dev', metavar='sd', required=True, dest='sd', help='standard deviation of library')
    parser.add_argument('-g', '--GTF', metavar='gtf', dest='gtf', required=False, default='blank',
                        help='reference gtf file to map reads to')
    parser.add_argument('--species', action="store", dest='species', default="mm9", required=False)
    parser.add_argument('--single-end', action='store_true', dest='single_end', default=False)
    parser.add_argument('--stranded', action='store_true', required=False, dest='strand', default=False, help="indicate if library contains strand information")
    parser.add_argument('--rmdup', action='store_true', required=False, dest='rmdup', default=False, help="remove duplicates")
    args = parser.parse_args()
    
    library_type = 'fr-unstranded'
    if args.strand: library_type = 'fr-secondstrand'
    
    tophat(args.date, args.fastq_files[0], args.single_end, \
           args.mean, args.sd, args.gtf, library_type, args.species)

if __name__ == "__main__":
    main(sys.argv)
