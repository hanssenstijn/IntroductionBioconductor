# Remove global environment
#-----------------------------------------------------#
rm(list = ls())

# Load tidyverse package
#-----------------------------------------------------#
library(Biostrings)
library(GenomicRanges)
library(IRanges)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("airway")
BiocManager::install("SingleCellExperiment")

library(BSgenome.Hsapiens.UCSC.hg38)
library(airway)
library(SingleCellExperiment)
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

# Hsapiens 
#-----------------------------------------------------#
chr14_range = GRanges("chr14", IRanges(1, seqlengths(Hsapiens)["chr14"]))
chr14_dna <- getSeq(Hsapiens, chr14_range)
letterFrequency(chr14_dna, "GC", as.prob=TRUE)

# Mouse CDS sequence
#-----------------------------------------------------#
url <- "ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz"
fl <- BiocFileCache::bfcrpath(rnames = url)
cds <- rtracklayer::import(fl, "fasta")
# Remove cds
pred1 <- width(cds) %% 3 == 0
table(pred1)
pred2 <- narrow(cds, 1, 3) == "ATG"
stops <- c("TAA", "TAG", "TGA")
pred3 <- narrow(cds, width(cds) - 2, width(cds)) %in% stops
table(pred1 & pred2 & pred3)
cds <- cds[ pred1 & pred2 & pred3 ]

# Subset & Visualize 
#-----------------------------------------------------#
hist(log10(width(cds)))
# Longest width
cds[ which.max(width(cds)) ]
# Check name
names(cds)[ which.max(width(cds)) ]
# Check GC content of each cds
gc <- letterFrequency(cds, "GC", as.prob=TRUE)
head(gc)
# Visualize
hist(gc)
plot( log10(width(cds)), gc, pch=".")

# Summarize codon usage
#-----------------------------------------------------#
AMINO_ACID_CODE
aa <- translate(cds)
codon_use <- letterFrequency(aa, names(AMINO_ACID_CODE))
head(codon_use)

# DNAStringSet
#-----------------------------------------------------#
mcols(cds) <- DataFrame(
    GC = gc[,"G|C"]
)
mcols(cds, use.names = FALSE)
mcols(cds[1:3], use.names = FALSE)

# Experiment
#-----------------------------------------------------#
data(airway)
airway

# Inspect data
#-----------------------------------------------------#
colData(airway)
airway[ , airway$dex == "untrt"]
# Assay() matrix-like container where rows represent featurs of interest
colSums(assay(airway))

# Mouse brain data set (rds)
#-----------------------------------------------------#
url <- "https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/manno_mouse.rds"
fl <- BiocFileCache::bfcrpath(rnames = url)
sce <- readRDS(fl)

# GenomicRanges
#-----------------------------------------------------#
