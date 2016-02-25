#!/bin/sh

set -e
#########################
# QN
../scripts/irap_qn --in ex1/pheno.tsv -o qn_test1.tsv -m ex1/sample2class.tsv
