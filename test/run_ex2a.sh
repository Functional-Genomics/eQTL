#!/bin/bash
# to override the corr_method run
# run_ex2.sh corr_method=panama
# or
# run_ex2.sh corr_method=peer
eqtl_pipeline conf=ex2.conf dna_rna_mapfile=ex1/map_file_test2.tsv  vcfs=vcf1/test.HG00096.vcf.gz\ test.HG00099.vcf.gz  $* 

