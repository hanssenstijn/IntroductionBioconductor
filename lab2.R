# Load packages
#-----------------------------------------------------#
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
library(Biostrings)
library(GenomicRanges)

# Inspect packages
#-----------------------------------------------------#
vignette(package="GenomicRanges")
vignette("GenomicRangesIntroduction")

# Create Genomic Ranges
#-----------------------------------------------------#
gr <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
    strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    score = 1:10,
    GC = seq(1, 0, length=10))
# Genomic coordinates located on the left hand side
# Metadata (annotation) right side
gr

# Inspect genomic coordinates
#-----------------------------------------------------#
# Genomic ranges
seqnames(gr)
ranges(gr)
strand(gr)
granges(gr)
# Annotations
mcols(gr)
mcols(gr)$score

# Splitting and combining GRanges objects
#-----------------------------------------------------#
sp <- split(gr, rep(1:2, each=5))
sp
# concatenate
c(sp[[1]], sp[[2]])
gr[2:3]

# Subset GRanges
#-----------------------------------------------------#
gr[2:3]
gr[2:3, "GC"]
singles <- split(gr, names(gr))
grMod <- gr
grMod[2] <- singles[[1]]
head(grMod, n=3)


# Basic interval operations for GRanges objects
#-----------------------------------------------------#
g <- gr[1:3]
g <- append(g, singles[[10]])
start(g)
flank(g, 10)

# Interval set operations for GRanges objects
#-----------------------------------------------------#
g2 <- head(gr, n=2)
union(g, g2)
setdiff(g, g2)

# Biostrings for DNA sequences
#-----------------------------------------------------#
data(phiX174Phage)
phiX174Phage
m <- consensusMatrix(phiX174Phage)[1:4,]
