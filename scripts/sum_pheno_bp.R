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
# You should have received a copy of the GNU General Public License if
# not, see <http://www.gnu.org/licenses/>.
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

usage <- "sum_pheno_bp.R -s summary_file -p phenotype -o output_tsv_file [-t title]"
option_list <- list(
  make_option(c("-o","--out"),type="character",default=NULL,dest="out",help="Output quantification file"),
  make_option(c("-s","--sum"),type="character",default=NULL,dest="sum_file",help="Quantification matrix"),
  make_option(c("-p","--pheno"),type="character",default=NULL,dest="pheno_file",help="Phenotype"),
  make_option(c("-t","--title"),type="character",default="",dest="title",help="Plot title"),
  make_option(c("--sig"),type="numeric",default=0.05,dest="sign.threshold",help=""))

multiple.options = list()
filenames <- c("sum_file","pheno_file") ;#filenames that must exist (if defined)

# check multiple options values
mandatory <- c("sum_file","out","pheno_file")

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
pheno <- fread(input=opt$pheno_file,sep = "\t" , header=TRUE,check.names=FALSE,data.table=FALSE,verbose=FALSE)
rownames(pheno) <- as.character(pheno[,1])


mean.pheno <- rowMeans(pheno[,-1,drop=FALSE])
print(head(mean.pheno))
pinfo("Data loaded.")
pinfo("Generating 2D plot ",opt$out,"...")

head(data)
egenes <- unique(as.character(data[signif,"geneID"]))
not.egenes <- unique(as.character(data[!signif,"geneID"]))

egenes.e <- mean.pheno[egenes]
not.egenes.e <- mean.pheno[not.egenes]
all.e <- mean.pheno[append(egenes,not.egenes)]

png(filename=opt$out,width=900,height=900,res=150)
par(bty="l")
boxplot(list("All genes"=all.e,"e-genes"=egenes.e,"Not e-genes"=not.egenes.e),outline=F,ylab="Expression",main=opt$title)
dev.off()
#
q(status=0)
save.image("l.Rdata")
load("l.Rdata")
