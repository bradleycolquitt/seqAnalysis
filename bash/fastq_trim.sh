#! /bin/bash

# $1 input fastq.gz
# $2 length of records in output

prefix=(${1//.fastq.gz/ })

gzip -dc $1 | awk -v outlength="$2" '
{
if (NR % 2 == 0) {
	    print substr($1, 1, outlength)
} else print 
}
' | gzip >  ${prefix[0]}_"$2"bp.fastq.gz
