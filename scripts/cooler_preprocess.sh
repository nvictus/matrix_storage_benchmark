#!/bin/bash

# Lexsort input pixels (assumes 3 column file: bin1, bin2, count)

INNAME=$1
OUTNAME=${INNAME/.tsv.gz/.3col.sorted.tsv.gz}

zcat $INNAME \
    | sort -k1,1n -k2,2n \
    | gzip -c \
    > $OUTNAME
