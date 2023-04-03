#!/bin/bash

less ${1}_detected.tsv | cut -f1 | tail -n +2 > names.txt

for val in $(cat names.txt)
do
	RScript plot_taxon.R $1 $val

done

rm names.txt
