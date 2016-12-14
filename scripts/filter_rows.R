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

rows.to.keep <- colnames(ref.df)

# 
cat("Reading tsv file...\n")
x <- read.table(tsv.file,sep="\t",header=T,check.names=F,nrow=1)
cat("Matrix with ",length(colnames(x))-1," data columns\n")

cat("Keeping  ",length(rows.to.keep)," rows\n")

nlines <- 0
stop = FALSE
# TODO: add checks
f = file(tsv.file, "r")
write.con = file(out.file, "w")
# Manual - print header
header = readLines(f, n = 1)
writeLines(header,write.con,useBytes=TRUE)

while(!stop) {
    next.line = unlist(strsplit(readLines(f, n = 1),split="\t",fixed=T),recursive=FALSE,use.names=F)
    
    ## Insert some if statement logic here
    if(length(next.line) == 0) {
        stop = TRUE
        close(f)
    } else {
        if ( next.line[1] %in% rows.to.keep ) {
            next.line <- paste(next.line,collapse="\t")
            writeLines(next.line,write.con,useBytes=TRUE)
            nlines <- nlines+1
            cat(".")
        }
    }
}
cat("\n")
#write.table(new.x,file=out.file,quote=F,sep="\t",row.names=F,col.names=T)
quit(status=0)

