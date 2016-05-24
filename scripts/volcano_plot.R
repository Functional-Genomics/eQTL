#!/bin/env Rscript
#; -*- mode: R;-*-
# =========================================================
# Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is based on code from iRAP (https://github.com/nunofonseca/irap)
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
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
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

# load a file with a quant. matrix
# returns NULL in case of failure
quant.load <- function(f,clean.cuff=FALSE) {
  tsv.data <- NULL

  # header is always present
  tsv.data <- qload.tsv(f,header=TRUE)
  if(is.null(tsv.data)) return(NULL);
  rownames(tsv.data) <- as.character(tsv.data[,1])
  tsv.data <- tsv.data[,-1,drop=FALSE]
  if (clean.cuff) {
    sel<-grep("^CUFF.*",rownames(tsv.data),perl=T,invert=TRUE)
    tsv.data <- tsv.data[sel,,drop=FALSE]
  }
  return(tsv.data)
}

# keep backwards compatibility by using read.table when data.table is not 
# available
qload.tsv <- function(f,header,comment.char="") {
  tsv.data <- NULL
  if (require("data.table",quietly=TRUE,character.only=TRUE) &&
      compareVersion(as.character(packageVersion("data.table")),"1.9.6")>=0) {
    library("data.table")
    if ( sum(grep(".gz$",f)) ) {
      f <- paste("zcat ",f,sep="")
    } else {
      f <- paste("cat ",f,sep="")
    }
    # not optimal, but faster than read.table
    if ( comment.char!="") {
      f <- paste(f," | grep -v \"^",comment.char,"\"",sep="")
    }
    tryCatch(tsv.data <- fread(input=f,sep = "\t", header=header,check.names=FALSE,data.table=FALSE),error=function(x) NULL)
  } else 
    tryCatch(tsv.data <- read.table(f,sep = "\t", header=header, comment.char=comment.char, quote = "\"",check.names=FALSE),error=function(x) NULL)
  return(tsv.data)
}

quantile_norm <- function(df,means=NULL){
  if ( ! is.data.frame(df) ) {
    perror("Expected a data frame")
  }  
  #
  if ( ! is.null(means)) {
    l1 <- nrow(df)
    l2 <- length(means)
    if ( l1 != l2 ) {
      pwarning("Number of rows in data frame (",l1,") does not match with the length of quantile normalized means vector (",l2,")")
      if ( l1 > l2 ) {
        perror("Unable to proceed")
      }
      # l1 <l2
      offset <- l2-l1
      means <- means[append(2,seq(offset+2,l2))]
      #pinfo(length(means),"==",l1)
    }  
  }
  print(dim(df))
  # increasing
  ranks <- apply(df,2,rank,ties.method="max")
  # sort: increasing
  if (is.null(means) ) {
    means <- apply(data.frame(apply(df, 2, sort)), 1, mean, na.rm=T)
  }
  df_qn<- apply(ranks, 2, quantile_norm_vect, means)
  rownames(df_qn) <- rownames(df)
  return(list(qn=df_qn,means=means))
}

quantile_norm_vect <- function(v,qn_values) {
  lv <- length(v)
  lqn <- length(qn_values)
  if ( lv!=lqn ) {
    perror("length of vector v (",lv,") is different from qn_means' length (",lqn,")")
  }
  
  p <- rank(v,ties.method="max")
  return(qn_values[p])
}

###############################################################

usage <- "volcano_plot -i summary_file -m mapping_tsv_file -o output_tsv_file"
option_list <- list(
  make_option(c("-o","--out"),type="character",default=NULL,help="Output quantification file"),
  make_option(c("-i","--in"),type="character",default=NULL,dest="tsv_file",help="Quantification matrix"),
  make_option(c("-t","--title"),type="character",default="",dest="title",help="Plot title"),
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

head(data$g_adj_pval)
plot(data$beta,-log10(data$g_adj_pval),
     main=opt$title,
     xlab="Beta",ylab="-log2(adj. pvalue)",pch=20,
     ylim=c(0,1+max(-log10(data$g_adj_pval))))
pos.sig.beta <- data$beta>0 & signif
neg.sig.beta <- data$beta<0 & signif
points(data$beta[pos.sig.beta],-log10(data$g_adj_pval[pos.sig.beta]),col="red",pch=20)
points(data$beta[neg.sig.beta],-log10(data$g_adj_pval[neg.sig.beta]),col="blue",pch=20)
abline(h=-log10(opt$sign.threshold),lty=2,col="grey")
abline(h=-log10(0.5),lty=2,col="black")
abline(v=0,lty=2,col="grey")
# add a few labels
sel.neg <- head(order(data$g_adj_pval[neg.sig.beta]),n=2)
sel.pos <- head(order(data$g_adj_pval[pos.sig.beta]),n=2)
#labels1 <- apply(data[neg.sig.beta,c("GenoName","GeneName")][sel.neg,],MARGIN=1,FUN=paste,collapse="/")
#labels2 <- apply(data[pos.sig.beta,c("GenoName","GeneName")][sel.pos,],MARGIN=1,FUN=paste,collapse="/")
par(xpd=TRUE)
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

q(status=0)
