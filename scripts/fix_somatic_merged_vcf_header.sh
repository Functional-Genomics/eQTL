#!/bin/bash

function usage()
{
        cat<<EOF

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


sample_name=$(echo $vcfin | cut -f 1 -d . ); subst_line=$(zcat $vcfin | grep '#CHROM' | sed "s/$/\t$sample_name/"); zcat $vcfin | sed "s/#CHROM.*$/$subst_line/g" |  bgzip -c > $vcfout.tmp && 

mv $vcfout.tmp $vcfout
rm -f $vcfout.tmp

exit 0



