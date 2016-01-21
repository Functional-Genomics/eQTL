#!/bin/bash

function usage()
{
        cat<<EOF

Usage:

fix_sanger_vcf_header.sh <filename.vcf.gz> <filename.fixedheader.vcf.gz>
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


sample_name=$(zcat $vcfin | grep '##SAMPLE=<ID=TUMOUR' | cut -f 8 -d '=' | cut -f 1 -d ','); zcat $vcfin | grep -v '##SAMPLE=<ID=NORMAL' | sed 's/TUMOUR/'$sample_name'/g' | cut -f 1,2,3,4,5,6,7,8,9,11 | bgzip -c > $vcfout.tmp &&

mv $vcfout.tmp $vcfout
rm -f $vcfout.tmp

exit 0



