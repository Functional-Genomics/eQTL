#!/bin/bash

function usage()
{
        cat<<EOF

This script just copy the input vcf file
Usage:

fix_merged_vcf_header.sh <filename.vcf.gz> <filename.fixedheader.vcf.gz>
EOF

}

vcfin=$1
vcfout=$2

if [ "$vcfin-" == "-" ]; then
    echo "ERROR: missing input file"
    usage
    exit 1
fi

if [ ! -e "$vcfin"  ]; then
    echo "ERROR: $vcfin file not found"
    exit 1
fi

if [ "$vcfout-" == "-" ]; then
    echo "ERROR: missing output file"
    usage
    exit 1
fi

cp $vcfin $vcfout

exit 0



