# Remove global environment
#-----------------------------------------------------#
rm(list = ls())

# Load tidyverse package
#-----------------------------------------------------#

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("tximport","tximportData","DESeq2","tximeta","pasilla"))

library("tximport")
library("tximportData")
library("readr")
library("DESeq2")
library("tximeta")
library("pasilla")

# Import data tximport
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

# Import data tximeta
#-----------------------------------------------------#
coldata <- samples
coldata$files <- files
coldata$names <- coldata$run
se <- tximeta(coldata)
ddsTxi <- DESeqDataSet(se, design = ~ condition)

# Count matrix input
#-----------------------------------------------------#
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
# Count matrix and column data are consistent
# It is absolutely critical that the columns of the count matrix and
# the rows of the column data (information about samples) are in the same order.
head(cts,2)
coldata
# Samen column order
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
# Add addtional feature data
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
