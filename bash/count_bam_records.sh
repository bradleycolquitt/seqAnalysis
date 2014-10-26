#! /bin/bash

samtools idxstats $1 | awk 'BEGIN{total=0}{total=$3+total}END{print total}'
