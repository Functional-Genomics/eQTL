#!/bin/bash

function usage() 
{
	cat<<EOF

Usage:

fix_vcf_header.sh <filename.vcf.gz> <filename.fixedheader.vcf> 
EOF
	
}

#no need to generate gz vcf in this step! TODO: change make file? """

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


zcat $vcfin | grep '##' | sed 's/ //g' > $vcfout &&
zcat $vcfin | grep -v '##' >> $vcfout


exit 0

