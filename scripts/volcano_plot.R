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

usage <- "volcano_plot -i summary_file -m mapping_tsv_file -o output_tsv_file"
option_list <- list(
  make_option(c("-o","--out"),type="character",default=NULL,help="Output quantification file"),
  make_option(c("-i","--in"),type="character",default=NULL,dest="tsv_file",help="Quantification matrix"),
  make_option(c("-t","--title"),type="character",default="",dest="title",help="Plot title"),
  make_option(c("--save_image"),action="store_true",dest="save.image",default=FALSE,help="Save the R image"),  
  make_option(c("-s","--sig"),type="numeric",default=0.05,dest="sign.threshold",help=""))

multiple.options = list()
filenames <- c("tsv_file") ;#filenames that must exist (if defined)

# check multiple options values
mandatory <- c("tsv_file","out")

#
args <- commandArgs(trailingOnly = TRUE)
opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory,args=args)

library(data.table)
pinfo("Loading ",opt$tsv_file)
data <- fread(input=opt$tsv_file,sep="\t",data.table=FALSE,header=T)
pinfo("Loading ",opt$tsv_file, " complete")
pinfo("#entries:",nrow(data))

if (nrow(data)==0) {
  cat("WARNING: No data to plot\n")
  png(opt$out,width=900,height=900,res=150)
  plot(0,0)
  text(0,0,"no data to plot")
  dev.off()
  q(status=0)
}

expected.cols <- c("g_adj_pval","beta")
if ( sum(! expected.cols %in% colnames(data))>0 ) {
  cat("ERROR: columns missing ",expected.cols[! expected.cols %in% colnames(data)],"\n")
  q(status=1)
}
signif <- data$g_adj_pval<=opt$sign.threshold

pinfo("Generating volcano plot...")
png(opt$out,width=900,height=900,res=150)
par(bty="l")

###
pinfo("pvalues=0: ",sum(data$g_adj_pval==0),"\n")
data$g_adj_pval[data$g_adj_pval==0] <- .Machine$double.xmin


plot(data$beta,-log10(data$g_adj_pval),
     main=opt$title,
     xlab="Beta",ylab="-log2(adj. pvalue)",pch=20,
     ylim=c(0,1+max(-log10(data$g_adj_pval))))

pos.sig.beta <- NULL
neg.sig.beta <- NULL
pos.sig.beta <- data$beta>0 & signif
neg.sig.beta <- data$beta<0 & signif
if ( sum(pos.sig.beta)>0 ) {
  points(data$beta[pos.sig.beta],-log10(data$g_adj_pval[pos.sig.beta]),col="red",pch=20)
}
if ( sum(neg.sig.beta)>0 ) {
  points(data$beta[neg.sig.beta],-log10(data$g_adj_pval[neg.sig.beta]),col="blue",pch=20)
}
abline(h=-log10(opt$sign.threshold),lty=2,col="grey")
abline(h=-log10(0.5),lty=2,col="black")
abline(v=0,lty=2,col="grey")
par(xpd=TRUE)
# add a few labels
if ( sum(neg.sig.beta)>0 ) {
  sel.neg <- head(order(data$g_adj_pval[neg.sig.beta]),n=2)
}
if ( sum(pos.sig.beta)>0 ) {
  sel.pos <- head(order(data$g_adj_pval[pos.sig.beta]),n=2)
}
#labels1 <- apply(data[neg.sig.beta,c("GenoName","GeneName")][sel.neg,],MARGIN=1,FUN=paste,collapse="/")
#labels2 <- apply(data[pos.sig.beta,c("GenoName","GeneName")][sel.pos,],MARGIN=1,FUN=paste,collapse="/")

text(x=max(data$beta),y=-log10(opt$sign.threshold),labels=paste("P<",opt$sign.threshold,sep=""),pos=1,cex=0.8)
text(x=max(data$beta),y=-log10(0.5),labels=paste("P<0.5",sep=""),pos=1,cex=0.8)
## text(data$beta[neg.sig.beta][sel.neg],-log10(data$g_adj_pval[neg.sig.beta][sel.neg]),
##      labels=labels1,
##      cex=0.5,
## pos=c(3,4))
## text(data$beta[pos.sig.beta][sel.pos],-log10(data$g_adj_pval[pos.sig.beta][sel.pos]),
##      labels=labels2,
##      cex=0.5,
##      pos=c(3,4))
dev.off()

if ( opt$save.image ) { 
    save.image(paste(opt$out,".Rdata",sep=""))
}
q(status=0)
