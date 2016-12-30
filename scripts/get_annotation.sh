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

awk '$3 == "$feature" {print $1,$4,$5,$10,$7}' $gtf | sed "s/chr//g;s/\;//g;s/ /\t/g;s/\"//g" > $out

exit 0

