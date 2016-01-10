#!/usr/bin/env Rscript
#; -*- mode: R;-*-

usage <- "genQuantBed.R annot_tsv_file quant_tsv_file out_filename"

args <- commandArgs(trailingOnly=TRUE)
if ( length(args) != 3 ) {
  cat("ERROR: Usage:",usage,"\n")
  quit(status=1);
}
annot.tsv.file <- args[1]
quant.tsv.file <- args[2]
out.file <- args[3]

library(data.table)

load.file <- function(f,header=TRUE) {
  cat("Loading ",f,"...\n")
  if ( ! file.exists(f) ) {
    cat("ERROR: File not found ",f,"\n")
    q(status=1)
  }
  x <- tryCatch(fread(input=f,header=header,"\t",
                      data.table=FALSE,showProgress=TRUE,colClasses="character"),error=function(x) return(NULL))
  if ( is.null(x) ) {
    cat("ERROR while loading ",f,"\n")
    q(status=1)
  }
  return(x)
}

annot <- load.file(annot.tsv.file,header=FALSE)
if ( ncol(annot)!=4 ) {
  cat("ERROR: expected 4 columns in ",annot.tsv.file,"\n")
  q(status=1)
}
colnames(annot) <- c("#Chr","start","end","ID")

quant <- load.file(quant.tsv.file,header=TRUE)
rownames(quant) <- quant[,1]
quant <- quant[,-1]
# sort quant based on the annot
miss <- (! rownames(quant) %in% as.character(annot$ID))
if ( sum(miss) >0 ) {
  cat("Inconsistent quant. and annot files - missing genes. E.g. ",head(annot$ID[miss]),"\n")
  q(status=1)
}
miss <- (! as.character(annot$ID) %in% rownames(quant))
if ( sum(miss) >0 ) {
  cat("Warning: Inconsistent quant. and annot files - missing genes in annotation. E.g. ",head(annot$ID[miss]),"\n")
}
sel <- as.character(annot$ID) %in% rownames(quant)
quant <- quant[as.character(annot$ID[sel]),]
nm <- cbind(annot[sel,],quant)
o <- order(as.character(nm[,1]),as.numeric(nm[,2]))
r <- write.table(nm[o,],file=out.file,col.names=TRUE,sep="\t",row.names=F,quote=F)
cat("Saved ",out.file,"\n")

q()


