#!/usr/bin/Rscript

# Select only the genes in the GTF file
args <- commandArgs(trailingOnly=TRUE)
tsv.file <- args[1]
gtf.file <- args[2]
out.file <- args[3]

library("data.table")
##################
# load the phenotype
cat("Loading expression file...\n")
ref.df<-fread(input=tsv.file,sep="\t",header=T,data.table=FALSE)
cat("Samples:",ncol(ref.df)-1,"\n")
cat("Entries:",nrow(ref.df),"\n")
if ( ncol(ref.df)<=1 ) {
  cat("ERROR: insufficient number of columns\n")
  q(status=1)
}
if ( colnames(ref.df)[1]=="" ) {
  colnames(ref.df)[1] <- "Gene"
}
s <- colnames(ref.df)[1]

rownames(ref.df) <- ref.df[,1]
# 

cat("Reading annot. file...\n")
x <- fread(input=gtf.file,sep="\t",header=FALSE,data.table=FALSE)
colnames(x) <- c("chr","start","end","id")
cat("Matrix with ",length(colnames(x))," data columns\n")

sel <- intersect(as.character(x$id),rownames(ref.df))
new.ref <- ref.df[sel,,drop=FALSE]
cat("New matrix with ",length(rownames(new.ref))," data columns\n")

write.table(new.ref,file=out.file,quote=F,sep="\t",row.names=F,col.names=T)
quit(status=0)
