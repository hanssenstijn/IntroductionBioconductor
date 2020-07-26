# Remove global environment
#-----------------------------------------------------#
rm(list = ls())

# Load tidyverse package
#-----------------------------------------------------#

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("tximport","tximportData",
                       "DESeq2","tximeta",
                       "pasilla","apeglm",
                       "BiocParallel"))

library("tximport")
library("tximportData")
library("readr")
library("DESeq2")
library("tximeta")
library("pasilla")
library("airway")
library("BiocParallel")
library("ggplot2")

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
# Same column order
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
# Add additional feature data
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# Htseq-count
#-----------------------------------------------------#
directory <- system.file("extdata", package="pasilla",
                         mustWork=TRUE)
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

# SummarizedExperiment input
#-----------------------------------------------------#
data("airway")
se <- airway
ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
ddsSE

# Pre-filtering
#-----------------------------------------------------#
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Factor levels
#-----------------------------------------------------#
# 2 different ways
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
# Second option: dds$condition <- relevel(dds$condition, ref = "untreated")
# Drop levels which are not being used anymore
dds$condition <- droplevels(dds$condition)

# Differential expression analysis
#-----------------------------------------------------#
dds <- DESeq(dds)
res <- results(dds)
res
# Specify
res <- results(dds, name="condition_treated_vs_untreated")
res <- results(dds, contrast=c("condition","treated","untreated"))

# Log fold change shrinkage for visualization and ranking
#-----------------------------------------------------#
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC

# Using parallelization
#-----------------------------------------------------#
register(SnowParam(4))

# p-values and adjusted p-values
#-----------------------------------------------------#
resOrdered <- res[order(res$pvalue),]
summary(res)
# How many adjusted p-values are less than 0.1
sum(res$padj < 0.1, na.rm=TRUE)
# Different cutoff value
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

# Exploring and exporting results
#-----------------------------------------------------#
# Shows the log2 fold changes attributable to a given variable over 
# the mean of normalized counts for all the samples in the DESeqDataSet
plotMA(res, ylim=c(-2,2))
# the shrunken log2 fold changes, which remove the noise associated 
# with log2 fold changes from low count genes
plotMA(resLFC, ylim=c(-2,2))
# After calling plotMA, one can use the function identify 
# to interactively detect the row number of individual genes 
# by clicking on the plot.
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# Plot counts
#-----------------------------------------------------#
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10(breaks=c(25,100,400))
mcols(res)$description

# Count data transformations
#-----------------------------------------------------#
# Extract the matrix of normalized values
vsd <- vst(dds, blind=FALSE)
