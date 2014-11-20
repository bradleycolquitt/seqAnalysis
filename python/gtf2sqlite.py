#! /usr/bin/env python

##############################################################################################
# Parse CompletedJobInfo.xml and GenerateFASTQRunStatistics.xml within Illumina run folder
# Export
#   - library ID <SampleID>
#   - run ID <Run>
#   - instrument <Instrument>
#   - cycle number <NumCycles>
#   - type of sequencing PE/SE (Calculate from number of children of <Read>)
#   - Total number of clusters <NumberOfClustersRaw>
#   - Passed number of clusters <NumberOfClustersPF>
##############################################################################################


import sys
import argparse
#import xml.etree.ElementTree as ET
#import MySQLdb as mdb
import sqlite3
import pdb
#import datetime
import traceback as tb

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='GTF file')
    args = parser.parse_args()

    gtf = open(args.gtf)

    ## Connect to db
    try:
        conn = sqlite3.connect("/media/data/db/anno.db", isolation_level=None)
        #pdb.set_trace()
        cur = conn.cursor()
        st = '''CREATE TABLE IF NOT EXISTS genes (build text,
                                                  chrom text,
                                                  element text,
                                                  start int,
                                                  end int,
                                                  strand text,
                                                  gene_id text,
                                                  transcript_id text
                                                  )'''
        cur.execute(st)
    except:
        #print "MySQLdb error %d: %s " % (e.args[0] + e.args[1])
        print tb.print_exc()


    chrom = "chr1"
    organism = "taeGut1"
    build = ""
    start = 0L
    end = 0L
    score = 0
    strand = "+"
    flag = ""
    element = "start_codon"
    gene_id = ""
    gene_id_discard = ""
    transcript_id = ""
    transcript_id_discard = ""

    record_parse = dict(chrom = "chr1",
                        build = "",
                        start = "",
                        end = "",
                        strand = "+",
                        element = "start_codon",
                        gene_id = "",
                        transcript_id = ""
                        )
    try:
        cur.execute("BEGIN TRANSACTION")
        for record in gtf:

            record_s = record.split()

            record_parse['chrom'] = "'" + record_s[0] + "'"
            record_parse['build'] = "'" + record_s[1] + "'"
            record_parse['element'] = "'" + record_s[2] + "'"
            record_parse['start'] = record_s[3]
            record_parse['end'] = record_s[4]
            record_parse['strand'] = "'" + record_s[6] + "'"
            record_parse['gene_id'] = "'" + record_s[9].split(";")[0] + "'"
            record_parse['transcript_id'] = "'" + record_s[11].split(";")[0] + "'"

            exec_string = "SELECT EXISTS(SELECT 1 FROM genes WHERE build=? AND transcript_id=?)"
            cur.execute(exec_string, (record_parse['build'], record_parse['transcript_id']))
            res = cur.fetchall()[0]
            if res[0] == 0:
                exec_string = "INSERT INTO genes (" + ", ".join(record_parse.keys()) + ") VALUES (" + ", ".join(record_parse.values()) + ")"
                cur.execute(exec_string)
            else:
                exec_string = "UPDATE genes SET " +\
                              ", ".join(["=".join(items) for items in record_parse.items()]) +\
                              " WHERE build=? AND transcript_id=?"
                cur.execute(exec_string, (record_parse['build'], record_parse['transcript_id']))
        conn.commit()
    except:
        conn.rollback()
        print tb.print_exc()
        #print "MySQLdb error %d: %s " % (e.args[0], e.args[1])
    #except:
    #    conn.rollback()
    #    print "MySQLdb programming error %d: %s " % (e.args[0], e.args[1])


if __name__ == "__main__":
    main(sys.argv)
