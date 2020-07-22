# Remove global environment
#-----------------------------------------------------#
rm(list = ls())

# Data packages
#-----------------------------------------------------#
"The name of an org package is always of the form org.<Sp>.<id>.db 
(e.g. org.Sc.sgd.db) where <Sp> is a 2-letter abbreviation of the organism
(e.g. Sc for Saccharomyces cerevisiae) and <id> is an abbreviation 
(in lower-case) describing the type of central identifier 
(e.g. sgd for gene identifiers assigned by the Saccharomyces Genome Database, or eg for Entrez gene ids)"

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(tidyverse)

# Check package
#-----------------------------------------------------#
keytypes(org.Hs.eg.db)
head(keys(org.Hs.eg.db, "SYMBOL"))

# Extract data from org.* package
#-----------------------------------------------------#
set.seed(123)
egid <- sample(keys(org.Hs.eg.db), 6)
mapIds(org.Hs.eg.db, egid, "SYMBOL", "ENTREZID")
# Or use select
AnnotationDbi::select(
    org.Hs.eg.db, egid, c("SYMBOL", "ENSEMBL", "GENENAME"), "ENTREZID"
)
egid <- "3812"
mapIds(org.Hs.eg.db, egid, "ENSEMBL", "ENTREZID")
# Many mapping between keys and columns
mapIds(
    org.Hs.eg.db, egid, "ENSEMBL", "ENTREZID",
    multiVals = "CharacterList"
)
AnnotationDbi::select(
    org.Hs.eg.db, egid, c("SYMBOL", "ENSEMBL"),
    multiVals = "CharacterList"
)
# Tidyverse results
egid <- keys(org.Hs.eg.db)    # all ENTREZIDs
mapIds(org.Hs.eg.db, egid, "SYMBOL", "ENTREZID") %>% 
    as_tibble() %>% 
    rownames_to_column("ENTREZID")
AnnotationDbi::select(
    org.Hs.eg.db, egid, c("SYMBOL", "GO", "GENENAME"), "ENTREZID"
) %>% as_tibble()
