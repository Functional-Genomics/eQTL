#!/usr/bin/Rscript


library("data.table")

args <- commandArgs(trailingOnly=TRUE)
tsv.file <- args[1]
tsv.ref.file <- args[2]
out.file <- args[3]

##################
# load the tsv file
cat("Loading reference file...\n")
ref.df<-read.table(tsv.ref.file,sep="\t",header=T,check.names=F,nrows=2)
ref.df <- ref.df[,-1]
cat("Samples:",ncol(ref.df),"\n")


# 
cat("Reading tsv file...\n")
x <- read.table(tsv.file,sep="\t",header=T,check.names=F,nrow=1)
cat("Matrix with ",length(colnames(x))-1," data columns\n")
if ( colnames(x)[1]=="" ) {
  colnames(x)[1] <- "Gene"
}
s <- colnames(x)[-1]
#head(s)
cols.to.keep <- append(TRUE,s %in% colnames(ref.df))
head(cols.to.keep)
new.x <- x[,append(colnames(x)[1],s[s %in% colnames(ref.df)])]
cat("New matrix with ",length(colnames(new.x))-1," data columns\n")

nlines <- 0
stop = FALSE
# TODO: add checks
f = file(tsv.file, "r")
write.con = file(out.file, "w")
while(!stop) {
    next.line = unlist(strsplit(readLines(f, n = 1),split="\t",fixed=T),recursive=FALSE,use.names=F)[cols.to.keep]
    
    ## Insert some if statement logic here
    if(length(next.line) == 0) {
        stop = TRUE
        close(f)
    } else {
        next.line <- paste(next.line,collapse="\t")
        writeLines(next.line,write.con,useBytes=TRUE)
        nlines <- nlines+1
        if ( nlines%%5000==0 ) {
            cat(".")
        }
    }
}
cat("\n")
#write.table(new.x,file=out.file,quote=F,sep="\t",row.names=F,col.names=T)
quit(status=0)

