#!/bin/bash

file=$1
column=$2
label2print=$3

# use bash and R to avoid further software dependencies
R --slave --quiet --vanilla <<EOF 
x<-read.table("$file",sep="\t",header=T,quote="")
cat("$label2print\t",sum(x[,$column]),"\n")
q(status=0)
EOF

exit 0
