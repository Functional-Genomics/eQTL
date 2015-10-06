#!/bin/bash

usage()
{
	USAGE=""" 
	get_annotation.sh <infile.gtf> <outfile.csv> """

	echo $USAGE
}

gtf=$1
csv=$2

awk '$3 == "gene" {print $1,$4,$5,$10}' $gtf | sed 's/chr//g' | sed 's/"//g' | sed 's/\;//g' | sed 's/ /\t/g' > $csv


