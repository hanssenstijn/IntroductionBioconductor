# Remove global environment
#-----------------------------------------------------#
rm(list = ls())

# Load tidyverse package
#-----------------------------------------------------#
library(Biostrings)
library(GenomicRanges)
library(IRanges)

# Import data
#-----------------------------------------------------#
data(phiX174Phage)
phiX174Phage

# Subset data
#-----------------------------------------------------#
m <- consensusMatrix(phiX174Phage)[1:4,]
polymorphic <- which(colSums(m != 0) > 1)
m[, polymorphic]

# Genomic ranges
#-----------------------------------------------------#
ir <- IRanges(start=c(10, 20, 30), width=5)
ir
# Identify adjacent ranges
flank(ir,3)
# GenmicRanges extends the notion of ranges to include
# features relevant to application of range in sequence analysis
gr <- GRanges(c("chr1", "chr1", "chr2"), ir, strand=c("+", "-", "+"))
gr

# Genomic ranges
#-----------------------------------------------------#
