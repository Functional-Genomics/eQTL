#!/bin/bash

function usage() 
{
	cat<<EOF

Usage:

fix_vcf_header.sh <filename.vcf.gz> <filename.fixedheader.vcf> 
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

named_pipe=$vcfout.tmp.pipe
#mkfifo $named_pipe
set -e
#to be optimised
zcat $vcfin | grep '^##' | sed 's/ //g' > $named_pipe && zcat $vcfin | grep -v '^##' >> $named_pipe && bgzip -c $named_pipe > $vcfout

#rm -f $named_pipe
exit 0

