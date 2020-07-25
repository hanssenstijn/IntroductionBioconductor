# Remove global environment
#-----------------------------------------------------#
rm(list = ls())

# Load tidyverse package
#-----------------------------------------------------#

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("tximport","tximportData","DESeq2"))

library("tximport")
library("tximportData")
library("readr")
library("DESeq2")

# Import data
#-----------------------------------------------------#
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$run
samples[,c("pop","center","run","condition")]

files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run
# Links transcripts to genes
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
# Quantification data
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
# Construct DESeqDataSet
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
