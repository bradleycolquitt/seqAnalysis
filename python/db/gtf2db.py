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
import xml.etree.ElementTree as ET
import MySQLdb as mdb
import pdb
import datetime

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='GTF file')
    args = parser.parse_args()

    gtf = open(args.gtf)

    ## Connect to db
    try:
        conn = mdb.connect('localhost', 'brad', 'Eu23ler1', 'annotations')
        #conn.autocommit(False)
        cur = conn.cursor()
    except mdb.Error, e:
        print "MySQLdb error %d: %s " % (e.args[0] + e.args[1])

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

            exec_string = "SELECT EXISTS(SELECT 1 FROM genes WHERE build={0} AND transcript_id={1})".format(record_parse['build'], record_parse['transcript_id'])
            cur.execute(exec_string)
            res = cur.fetchall()[0]
            if res[0] == 0:
                exec_string = "INSERT INTO genes (" + ", ".join(record_parse.keys()) + ") VALUES (" + ", ".join(record_parse.values()) + ")"
                cur.execute(exec_string)
            else:
                exec_string = "UPDATE genes SET " +\
                              ", ".join(["=".join(items) for items in record_parse.items()]) +\
                              " WHERE build={0} AND transcript_id={1}".format(record_parse['build'], record_parse['transcript_id'])
                cur.execute(exec_string)
        conn.commit()
    except mdb.OperationalError, e:
        conn.rollback()
        print "MySQLdb error %d: %s " % (e.args[0], e.args[1])
    except mdb.ProgrammingError, e:
        conn.rollback()
        print "MySQLdb programming error %d: %s " % (e.args[0], e.args[1])


if __name__ == "__main__":
    main(sys.argv)
