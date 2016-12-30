#!/bin/bash

function usage() {

    cat<<EOF    
Usage: get_annotation.sh <infile.gtf> feature <outfile.csv>
EOF

}

gtf=$1
feature=$2
out=$3

if [ "$gtf-" == "-" ]; then
    echo "ERROR: missing infile.gtf"
    usage
    exit 1
fi
if [ ! -e "$gtf"  ]; then
    echo "ERROR: $gtf file not found"
    exit 1
fi
if [ "$out-" == "-" ]; then
    echo "ERROR: missing outfile"
    usage
    exit 1
fi

if [ "-$feature" == "-gene" ]; then
    awk '$3 == "gene" {print $1,$4,$5,$10,$7}' $gtf | sed "s/chr//g;s/\;//g;s/ /\t/g;s/\"//g" > $out
else
    if [ "-$feature" == "-transcript" ]; then
	awk '$3 == "transcript" {print $1,$4,$5,$12,$7}' $gtf | sed "s/chr//g;s/\;//g;s/ /\t/g;s/\"//g" > $out
    else
	echo "INVALID feature $feature: supported features - gene, transcript" > /dev/stderr
    fi
fi


exit 0

