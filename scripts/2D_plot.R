#!/bin/env Rscript
#; -*- mode: R;-*-
# =========================================================
# Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# if not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================

if (!require("optparse")) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("optparse")
}
suppressPackageStartupMessages(library("optparse"))

###############################################################
# Functions from iRAP (copy to avoid an extra dependency)
perror <- function(...) {
  cat(paste("[ERROR] ",...,"\n",sep=""),file=stderr())
}
pinfo <- function(...) {
  cat(paste("[INFO] ",...,"\n",sep=""),file=stdout())
}

myParseArgs <- function(usage,option_list,filenames.exist=NULL,multiple.options=NULL,mandatory=NULL,...) {

  # get command line options, if help option encountered print help and exit,
  # otherwise if options not found on command line then set defaults,
  parser <- OptionParser(usage = usage, option_list=option_list)
  opt <- parse_args(parser,...)
  
  for ( m in mandatory ) {
    if ( is.null(opt[[m]]) ) {
        perror("Parameter ",m," needs to be defined")
        q(status=1)
    }
  }  
  for ( p in filenames.exist ) {
    if (! is.null(opt[[p]]) ) {
      if (! file.exists(opt[[p]]) ) {
        perror("File ",opt[[p]]," not found")
        q(status=1)
      }
    }
  }

  for ( op in names(multiple.options) ) {
    if ( ! opt[[op]] %in% multiple.options[[op]] ) {
      perror("Invalid value ",opt[[op]]," for option ",op)
      q(status=1)
    }
  }
  return(opt)
}



###############################################################

usage <- "2D_plot.R -s summary_file -p phenotype_info -c chr_sizes -o output_tsv_file [-t title -s sign. threshold]"
option_list <- list(
  make_option(c("-o","--out"),type="character",default="",dest="out",help="Output image"),
  make_option(c("-s","--sum"),type="character",default=NULL,dest="sum_file",help="Summary file"),
  make_option(c("-p","--pheno"),type="character",default=NULL,dest="pheno_file",help="Phenotype information"),
  make_option(c("-c","--chr"),type="character",default=NULL,dest="chr_file",help="Chr. sizes"),
  make_option(c("-x","--xlab"),type="character",default="",dest="xlab",help="x-label"),
  make_option(c("-y","--ylab"),type="character",default="",dest="ylab",help="y-label"),
  make_option(c("-t","--title"),type="character",default="",dest="title",help="Plot title"),
  make_option(c("--save_image"),action="store_true",dest="save.image",default=FALSE,help="Save the R image"),  
  make_option(c("--sig"),type="numeric",default=0.05,dest="sign.threshold",help=""))

multiple.options = list()
filenames <- c("sum_file","pheno_file","chr_file") ;#filenames that must exist (if defined)

# check multiple options values
mandatory <- c("sum_file","out","pheno_file","chr_file")

#
args <- commandArgs(trailingOnly = TRUE)
opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory,args=args)

library(data.table)


makeTransparent<-function(someColor, alpha=180) {
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#opt <- list()
#opt$sum_file <- "test_files/summary.tsv"
#opt$chr_file <- "ex1/chr_sizes.txt"
#opt$pheno_file <- "test_files/genes_pos.tsv"

pinfo("Loading ",opt$tsv_file)
data <- fread(input=opt$sum_file,sep="\t",data.table=FALSE,header=T)
pinfo("Loading ",opt$sum_file, " complete")
pinfo("#entries:",nrow(data))

if (nrow(data)==0) {
  cat("WARNING: No data to plot...exiting.\n")
  png(filename=opt$out,width=900,height=900,res=150)
  plot(0,0)
  text(0,0,"no data to plot")
  dev.off()
  q(status=0)
}

expected.cols <- c("g_adj_pval","beta","geneID","chrom","pos")
if ( sum(! expected.cols %in% colnames(data))>0 ) {
  cat("ERROR: columns missing ",expected.cols[! expected.cols %in% colnames(data)],"\n")
  q(status=1)
}
signif <- data$g_adj_pval<=opt$sign.threshold

#
pheno.annot <- fread(input=opt$pheno_file,sep = "\t" , header=FALSE,check.names=FALSE,data.table=FALSE,verbose=FALSE)
colnames(pheno.annot) <- c("Gene Chr","Gene Start","Gene End","Name")
#colnames(pheno.annot) <- c("Gene Chr","Gene Start","Gene End","Name","Strand")
rownames(pheno.annot) <- pheno.annot[,"Name"]
head(pheno.annot)

#
chr.length <- fread(input=opt$chr_file,sep = "\t", header=FALSE,check.names=FALSE,data.table=FALSE,verbose=FALSE)
rownames(chr.length) <- as.character(chr.length[,1])
colnames(chr.length) <- c("chr","size")
# chr of interest are the ones appearing in the pheno file
sel.chr <- intersect(rownames(chr.length),as.character(pheno.annot$"Gene Chr"))
chr.length <- chr.length[sel.chr,]

sel2 <- as.character(data$chr) %in% sel.chr
if ( sum(sel2) != nrow(data) ) {
  cat("WARNING: chromosomes in the summary file not found in the phenotype annotation file\n")   
}
data <- data[sel2,,drop=FALSE]

# Add the gene loci to the summary table
not.found <- as.character(data$geneID)[!as.character(data$geneID) %in% rownames(pheno.annot)]
not.found
if ( length(not.found) > 0 ) {
  cat("ERROR: genes in the summary file not in ",opt$pheno.file,"...e.g.,",head(not.found),"\n")
  q(status=1)
}
data <- cbind(data,pheno.annot[as.character(data$geneID),,drop=FALSE])
#head(data)
pinfo("Data loaded.")
pinfo("Generating 2D plot ",opt$out,"...")

#
n <- rownames(chr.length)
chr.length <- as.numeric(chr.length[n,2])
names(chr.length) <- n
chr.offset <- cumsum(as.numeric(chr.length))-chr.length
names(chr.offset) <- names(chr.length)


# exclude MT?...no

# Genotype names
#loci <- apply(geno.annot[,c(1,2)],c(1,2),as.character)
#loci <- gsub("\\s+","",loci)
#loci <- apply(loci,1,paste,collapse="-")
#loci2Name <- as.character(geno.annot$Name)
#names(loci2Name) <- loci
#head(loci2Name)

sum.loci <- apply(data[,c("chrom","pos")],c(1,2),as.character)
sum.loci <- gsub("\\s+","",sum.loci)
sum.loci <- apply(sum.loci,1,paste,collapse="-")
#head(sum.loci)

#fqtl.sum$GenoName <- loci2Name[fqtl.sum.loci]
#head(fqtl.sum)
geno.pos <- chr.offset[as.character(data$chrom)]+ data$pos

#pheno.pos <- geno.pos
pheno.pos <- chr.offset[data$"Gene Chr"]+data$"Gene Start"
col <- rep("black",length(pheno.pos))
col[data$beta<0] <- "blue"
col[data$beta>0] <- "red"
#par(
#length(geno.pos)
#length(pheno.pos)

png(filename=opt$out,width=900,height=900,res=150)
#par(bty="l",mar=c(4,4,4,4),mgp=c(3,0.5,0.5))
par(bty="l",mar=c(4,4,4,4),mgp=c(2.5,0.2,0))
plot(x=geno.pos,y=pheno.pos,
     xlab=opt$xlab,
     ylab=opt$ylab,
     main=opt$title,
     col=col,pch=20,
     xlim=c(0,max(chr.offset)+chr.length[1]),
     ylim=c(0,max(chr.offset)+chr.length[1]),
     axes=F,yaxt="n",xaxt="n",
     asp=1,xaxs="i",yaxs="i")
#chr.offset
#head(cumsum(chr.length))
rects <- data.frame(xl=cumsum(chr.length)-chr.length,
                    bt=rep(0,length(chr.offset)),
                    xr=cumsum(chr.length),
                    xt=sum(chr.length))
rownames(rects) <- paste("Chr",names(chr.offset),sep="")
rects.col <- c(makeTransparent("white"),makeTransparent("lightgrey"))
#rects

for (r in seq(1,nrow(rects))) {
  rect(rects[r,1],rects[r,2],rects[r,3],rects[r,4],col=rects.col[r%%2+1],border="transparent")
}
for (r in seq(1,nrow(rects))) {
  segments(0,rects[r,1],sum(chr.length),rects[r,1],col="grey")
}
segments(0,sum(chr.length),sum(chr.length),sum(chr.length),col="grey")
segments(sum(chr.length),0,sum(chr.length),sum(chr.length),col="grey")

# replot
points(x=geno.pos,y=pheno.pos,
     col=col,pch=20)
par(xpd=TRUE)
axis(1, at=append(0,cumsum(chr.length)), labels=FALSE,col="darkgrey",
          asp=1,xaxs="i",yaxs="i")
axis(1, at=cumsum(chr.length)-chr.length/2, labels=rownames(rects), las=2,col.axis="darkgrey",tck=0,col="darkgrey",cex.axis=0.6)

axis(2, at=append(0,cumsum(chr.length)), labels=FALSE,col="darkgrey",
     asp=1,xaxs="i",yaxs="i")
axis(2, at=cumsum(chr.length)-chr.length/2, labels=rownames(rects), las=2,col="darkgrey",tck=0,col.axis="darkgrey",cex.axis=0.6)

#legend("right",inset=c(-0.1,0),title="Effect\nsize",c(">0","<0","0"),pch=20,col=c("red","blue","black"),bty="n")
legend("right",inset=c(-0.1,0),title="Effect\nsize",c(">0","<0"),pch=20,col=c("red","blue"),bty="n")
dev.off()

if ( opt$save.image ) { 
    save.image(paste(opt$out,".Rdata",sep=""))
}
q(status=0)
