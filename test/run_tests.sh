#!/bin/sh

set -e
mkdir -p test_out
#########################
# QN
../scripts/irap_qn --in ex1/pheno.tsv -o test_out/qn_test1.tsv -m ex1/sample2class.tsv

../scripts/irap_qn --in files/expr1.fpkm.tsv.gz -o test_out/qn_expr1  -m files/expr1_map2class.tsv
