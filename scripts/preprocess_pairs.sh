#!/bin/bash

COOLER_PATH=''

time $COOLER_PATH/scripts/sort_contacts_txt.py \
	--chrom1 3 --pos1 4 --strand1 2 \
	--chrom2 7 --pos2 8 --strand2 6 \
	--out GM12878-MboI-contacts.txt.gz \
	b37.chrom.sizes $1
