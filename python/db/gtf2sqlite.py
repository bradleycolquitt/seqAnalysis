#! /usr/bin/env python

##############################################################################################
# Insert GTF fields into annotation DB
##############################################################################################


import sys
import argparse
import sqlite3
import pdb
import datetime
import re
import traceback as tb

class AnnoDB:
    def __init__(self, gtf):
        self.gtf_name = gtf
        self.gtf = open(self.gtf_name)
        self.build = ""
        now = datetime.datetime.now()
        self.date = "%s-%s-%s" % (now.year, now.month, now.day)
        self.nrows = 0
        self.conn = sqlite3.connect("/media/data/db/anno2.db", isolation_level=None)
        self.cur = self.conn.cursor()

        script = '''CREATE TABLE IF NOT EXISTS genes (build text,
                                                  chrom text,
                                                  element text,
                                                  start int,
                                                  end int,
                                                  strand text,
                                                  gene_id text,
                                                  transcript_id text,
                                                  gene_name text
                                                  );
                    CREATE TABLE IF NOT EXISTS builds (build text,
                                                  date_entered text,
                                                  nrows int
                                                  );
                 '''

        self.cur.executescript(script)

    def insert_anno(self):
        record_parse = dict(chrom = "chr1",
                            build = "",
                            start = "",
                            end = "",
                            strand = 0,
                            element = "start_codon",
                            gene_id = "",
                            transcript_id = "",
                            gene_name = ""
                            )

        try:
            self.cur.execute("BEGIN TRANSACTION")
            for record in self.gtf:
                self.nrows += 1
                if self.nrows % 10000 == 0:
                    print "processed: " + str(self.nrows)
                    self.cur.execute("COMMIT")
                    self.cur.execute("BEGIN TRANSACTION")
                record_s = record.split()
                record_parse['chrom'] = "'" + record_s[0] + "'"
                record_parse['build'] = "'" + record_s[1] + "'"
                record_parse['element'] = "'" + record_s[2] + "'"
                record_parse['start'] = record_s[3]
                record_parse['end'] = record_s[4]
                record_parse['strand'] = "0" if record_s[6] == "+" else "1"
                record_parse['gene_id'] = "'" + record_s[9].split(";")[0] + "'"
                record_parse['transcript_id'] = "'" + record_s[11].split(";")[0] + "'"

                if re.search("gene_name", record):
                    record_parse['gene_name'] = "'" + record_s[13].split(";")[0] + "'"

                # Check if record exists
                exec_string = '''SELECT EXISTS(SELECT 1 FROM genes
                                               WHERE build=? AND transcript_id=?)'''
                self.cur.execute(exec_string, (record_parse['build'],
                                               record_parse['transcript_id']))
                res = self.cur.fetchall()[0]
                if res[0] == 0:
                    # If record doesn't exist, insert into db
                    exec_string = "INSERT INTO genes (" +\
                                  ", ".join(record_parse.keys()) +\
                                  ") VALUES (" + ", ".join(record_parse.values()) + ")"
                    self.cur.execute(exec_string)
                else:
                    # Else update with current info
                    exec_string = "UPDATE genes SET " +\
                                  ", ".join(["=".join(items) for items in record_parse.items()]) +\
                                  " WHERE build=? AND transcript_id=?"
                    self.cur.execute(exec_string, (record_parse['build'],
                                                   record_parse['transcript_id']))
            self.conn.commit()
            self.build = record_parse['build']
        except:
            print tb.print_exc()
            self.conn.rollback()

    def update_build_table(self):
        try:
            exec_string = '''SELECT EXISTS(SELECT 1 FROM builds
                                               WHERE build=? AND date_entered=?)'''
            self.cur.execute(exec_string, (self.build, self.date))
            res = self.cur.fetchall()[0]
            if res[0] == 0:
                self.cur.execute('''INSERT INTO builds VALUES (?, ?, ?)''', (self.build,
                                                                             self.date,
                                                                             self.nrows))
            else:
                self.cur.execute('''UPDATE builds SET build=?, date_entered=?, nrows=?''', (self.build,
                                                                                   self.date,
                                                                                   self.nrows))
            self.conn.commit()
        except:
            self.conn.rollback()
            print tb.print_exc()

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='GTF file')
    args = parser.parse_args()

    annodb = AnnoDB(args.gtf)

    print "Inserting GTF records..."
    annodb.insert_anno()
    print "Updating build table..."
    annodb.update_build_table()



if __name__ == "__main__":
    main(sys.argv)
