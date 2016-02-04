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
# TODO: check if bed has 5 columns
col=`head -n 1 $bed_file|wc -w`
if [ $col -lt 5 ]; then
    echo "ERROR:bed file should have 5 columns" > /dev/stderr
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
set -e
tmp_file=`mktemp`
tmp_file2=`mktemp`
tmp_file3=`mktemp`


cut -f 1,2,3 $bed_file > $tmp_file3
cut -f 1,2,3,4,5 $bed_file > $tmp_file
for vcf in $vcfs; do
    bedtools intersect -c -a $tmp_file3 -b $vcf | cut -f 4- |paste $tmp_file - > $tmp_file2
    mv $tmp_file2 $tmp_file
done
echo chr start end name score $vcfs | tr " " "\t"
cat $tmp_file
rm -f $tmp_file $tmp_file3 $tmp_file2
exit 0
