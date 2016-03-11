#!/usr/bin/Rscript

version <- "V1"

args <- commandArgs(trailingOnly=TRUE)
tsv.file <- args[1]
tsv.ref.file <- args[2]
out.file <- args[3]

#tsv.file <- "data/freeze3_v2.tophat2.lib.gene.fpkm.gl.donor.tsv.gz"
#tsv.file2 <- "fqtl_fused_genes.tsv"
##################
# load the fusions
cat("Loading reference file...\n")
ref.df<-read.table(tsv.ref.file,sep="\t",header=T,check.names=F,nrows=2)
ref.df <- ref.df[,-1]
cat("Samples:",ncol(ref.df),"\n")

# 
cat("Reading tsv file...\n")
x <- read.table(tsv.file,sep="\t",header=T,check.names=F)
cat("Matrix with ",length(colnames(x))-1," data columns\n")
if ( colnames(x)[1]=="" ) {
  colnames(x)[1] <- "Gene"
}
s <- colnames(x)[-1]
head(s)
new.x <- x[,append(colnames(x)[1],s[s %in% colnames(ref.df)])]
cat("New matrix with ",length(colnames(new.x))-1," data columns\n")

write.table(new.x,file=out.file,quote=F,sep="\t",row.names=F,col.names=T)
quit(status=0)
