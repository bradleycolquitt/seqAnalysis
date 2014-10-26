#! /usr/bin/env python 

##############################################################################################
# Load metadata for aligment file to sample_db
# Assumes the following directory structure:
#   
#     project/aligner_index/libN/seqN/bamfile.bam
#   
# idsequencing should be determined by date and libN
# Export
#   - idsequencing - from directory
#   - filename - provided as argument
#   - aligner_index - from directory
#   - align_datetime - from file timestamp
#   - aligner - from BAM header
#   - total_reads - read from 'idxstats'

##############################################################################################


import sys
import os.path
import datetime
import MySQLdb as mdb
import pdb
import datetime
import pysam

def main(argv):

    fullname = os.path.abspath(argv[1])
    bamfile = pysam.Samfile(fullname)
    header = bamfile.header
    fullname_s = fullname.split("/")

    # Populate info dict
    info = dict()    
    info['idsequencing'] = fullname_s[-2].split("seq")[1]
    info['filename'] = os.path.basename(fullname)
    info['aligner_index'] = fullname_s[-4]
    bam_datetime = b = datetime.datetime.fromtimestamp(os.stat(fullname).st_mtime)
    info['align_datetime'] = bam_datetime.strftime("%Y-%m-%d %H:%M:%S")
    info['aligner'] = header['PG'][0]['PN']
    info['command'] = header['PG'][0]['cl']

    # Compute total number of aligned reads
    stats = pysam.idxstats(fullname)
    stats = [el.split("\t") for el in stats]
    total_reads = 0
    for el in stats:
        total_reads += int(el[2])
    info['total_reads'] = total_reads
 
    # Format dict entries for MySQL
    for i in info.iterkeys():
        info[i] = "'" + str(info[i]) + "'"    

    ## Connect to db
    try:
        conn = mdb.connect('localhost', 'brad', 'Eu23ler1', 'sample_db')
        cur = conn.cursor()
    except mdb.Error, e:
        print "MySQLdb error %d: %s " % (e.args[0] + e.args[1])

    #Write to db            
    ## Check if record exists using idsequencing, aligner_index, and aligner
    try:
        #pdb.set_trace()
        # Get idsequencing from idlibrary and date
        exec_string = "SELECT EXISTS(SELECT 1 FROM alignments WHERE align_datetime={0})".format(info['align_datetime'])
        cur.execute(exec_string)
        res = cur.fetchall()[0]
        if res[0] == 0:
            exec_string = "INSERT INTO alignments (" + ", ".join(info.keys()) + ") VALUES (" + ", ".join(info.values()) + ")"
            cur.execute(exec_string)
        else:
            exec_string = "UPDATE alignments SET " +\
                              ", ".join(["=".join(items) for items in info.items()] + ["=".join(items) for items in info.items()]) +\
                              " WHERE idsequencing={0} AND align_datetime={1}".format(info['idsequencing'], info['align_datetime'])
            cur.execute(exec_string)
        conn.commit()
    except mdb.OperationalError, e:
        print e
        
    #Close connection
    conn.close()

if __name__ == "__main__":
    main(sys.argv)
