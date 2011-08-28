import os, sys
import pysam
import bam2bed
import sam
import argparse
from subprocess import Popen

def tophat(inputs, output, mean, sd, gtf, library_type):
    if os.path.exists(output):
        dec = raw_input("Output file exists. Bypass mapping [y/n]? ")
        if dec == 'y': return
    if not os.path.exists(output): 
        os.mkdir(output)
    if gtf != 'blank': 
        cmd_args = ['tophat', '-p', '4', '-o', output, 
                    '-r', mean, '--mate-std-dev', sd, 
                    '-G', gtf, '--library-type', library_type, 
                    'mm9', inputs[0], inputs[1]]
    else:
        cmd_args = ['tophat', '-p', '6', '-o', output, 
                    '-r', mean, '--mate-std-dev', sd,
                    '--library-type', library_type,
                    'mm9', inputs[0], inputs[1]]
    print "Mapping with tophat: " + " ".join(cmd_args[1:])
    tophat = Popen(cmd_args)
    tophat.wait()

def main(argv):
    bam_dir = "/media/storage2/data/rna/tophat/"
    bed_dir = "/media/storage2/data/rna/bed/" 
    parser = argparse.ArgumentParser(description="Run tophat on given fastq files.\nOutputs to " + bam_dir)
    parser.add_argument('-i', '--input_fastq', nargs=2, metavar='fastq', required=True, 
                        dest='fastq_files', help='read1 and read3 of library')
    parser.add_argument('-n', '--index', dest='index', required=True, help='index number of library')
    parser.add_argument('-r', '--mate-inner-dist', metavar='mean', required=True, 
                        dest='mean', help='mean size of library (minus adaptors)')
    parser.add_argument('-s', '--mate-std-dev', metavar='sd', required=True, dest='sd', help='standard deviation of library')
    parser.add_argument('-g', '--GTF', metavar='gtf', dest='gtf', default='blank',
                        help='reference gtf file to map reads to')
    parser.add_argument('--stranded', action='store_true', required=False, dest='strand', default=False, help="indicate if library contains strand information")
    args = parser.parse_args()
     
    input_prefix = args.fastq_files[0].split('_')[0] + "_" + args.index
    output_path =  bam_dir + input_prefix 
    library_type = 'fr-unstranded'
    if args.strand: library_type = 'fr-secondstrand'
    tophat(args.fastq_files, output_path, args.mean, args.sd, args.gtf, library_type)
    
    #bamfile = bam_dir + "/accepted_hits.bam"
    bamfile = bam_dir + input_prefix + "/accepted_hits.bam"
    #sam.sam2bam(samfile, bamfile) 
    sam.proc(bamfile)
    
    #bam = output_path + input_prefix + "a.bam"
    #bed = bed_dir + input_prefix + ".bed"
    #sam.bam2bed(bam, bed)

    
if __name__ == "__main__":
    main(sys.argv)
