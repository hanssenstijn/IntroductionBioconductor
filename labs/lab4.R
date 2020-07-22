# Remove global environment
#-----------------------------------------------------#
rm(list = ls())

# Load tidyverse package
#-----------------------------------------------------#
library(Biostrings)
library(GenomicRanges)
library(IRanges)
library(Rsamtools)
library(GenomicAlignments)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("airway")
BiocManager::install("SingleCellExperiment")
BiocManager::install("RNAseqData.HNRNPC.bam.chr14")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("VariantAnnotation")

library(BSgenome.Hsapiens.UCSC.hg38)
library(airway)
library(SingleCellExperiment)
library('RNAseqData.HNRNPC.bam.chr14')
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(RNAseqData.HNRNPC.bam.chr14)
library(BiocParallel)

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
gr <- GRanges(c("chr1:10-14:+", "chr1:20-24:+", "chr1:22-26:+"))
shift(gr, 1) 
range(gr) 
reduce(gr)   
coverage(gr)
setdiff(range(gr), gr)    
# SNP overlap genes
genes <- GRanges(c("chr1:30-40:+", "chr1:60-70:-"))
snps <- GRanges(c("chr1:35", "chr1:60", "chr1:45"))
countOverlaps(snps, genes) > 0
# Gene nearest to the regulatory region
reg <- GRanges(c("chr1:50-55", "chr1:75-80"))
nearest(reg, genes)
precede(reg, genes)

# Aligned reads
#-----------------------------------------------------#
roi <- GRanges("chr14", IRanges(19653773, width=1)) 
bf <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[[1]], asMates=TRUE)
paln <- readGAlignmentsList(bf)
j <- summarizeJunctions(paln, with.revmap=TRUE)
j_overlap <- j[j %over% roi]
paln[j_overlap$revmap[[1]]]

# Variants
#-----------------------------------------------------#
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")
seqlevels(vcf) <- "chr22"
coding <- locateVariants(rowRanges(vcf),
                         TxDb.Hsapiens.UCSC.hg19.knownGene,
                         CodingVariants())
head(coding)

# Count the number of reads overlapping exons grouped into genes
#-----------------------------------------------------#
# Identify the regions of interest
exByGn <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
## only chromosome 14
seqlevels(exByGn, pruning.mode="coarse") = "chr14"
# sample BAM files
length(RNAseqData.HNRNPC.bam.chr14_BAMFILES)
# Summarize overlaps
olaps <- summarizeOverlaps(exByGn, RNAseqData.HNRNPC.bam.chr14_BAMFILES)
olaps
head(assay(olaps))
colSums(assay(olaps))  
plot(sum(width(olaps)), rowMeans(assay(olaps)), log="xy")

# Session info
#-----------------------------------------------------#
sessionInfo()