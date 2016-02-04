#!/bin/bash
# usage: vcfs2matrix.sh bed.file vcf1 vcf2 ...
bed_file=$1
usage="vcfs2matrix.sh bed.file [vcf1 vcf2 ...]"
if [ "x$bed_file" == "x" ] ; then
    echo "ERROR: $usage" > /dev/stderr
    exit 1
fi
if [ ! -e $bed_file ] ; then
    echo "ERROR: file $bed_file not found." > /dev/stderr
    echo "ERROR: $usage" > /dev/stderr
    exit 1
fi
shift 1
vcfs=$*
if [ "x$vcfs" == "x" ]; then
    # read the vcfs from stdin
    read -n vcfs
fi
#echo VCFS=$vcfs > /stderr
if [ "x$vcfs" == "x" ]; then
    echo "VCFs not provided" > /dev/stderr
    echo "ERROR: $usage" > /dev/stderr
    exit 1
fi
tmp_file=`mktemp`
tmp_file2=`mktemp`

set -e
cut -f 1,2,3 $bed_file > $tmp_file
col=`head -n 1 $bed_file|wc -w`
for vcf in $vcfs; do
    cut -f 1,2,3 $bed_file | bedtools intersect -c -a stdin -b $vcf | cut -f $col |paste $tmp_file - > $tmp_file2
    mv $tmp_file2 $tmp_file
done
echo chr start end $vcfs | tr " " "\t"
cat $tmp_file
rm -f $tmp_file
exit 0
