#!/bin/bash

function usage() {

    cat<<EOF    
Usage: get_annotation.sh <infile.gtf> <outfile.csv>
EOF

}

gtf=$1
out=$2

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

awk '$3 == "gene" {print $1,$4,$5,$10}' $gtf | sed "s/chr//g;s/\;//g;s/ /\t/g" > $out

exit 0

