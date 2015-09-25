#!/bin/bash


# eqtl

#
# Peer
#


#
# install limix
#

limix_path=
#
pythonfiles2fix=scripts/geno_preprocessing.py
for f in $pythonfiles2fix; do
    sed -i "s/^limix_path=.*/limix_path=$limix_path/" $f
done
