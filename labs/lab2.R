# Load tidyverse package
#-----------------------------------------------------#
library(readr)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)

# Import data
#-----------------------------------------------------#
fname  <- file.choose()
stopifnot(file.exists(fname))
brfss <- read_csv(fname)
brfss

# Using tidyverse
#-----------------------------------------------------#
# Operations in the tidyverse are ‘piped’ using the %>% symbol
# from one representation to another
brfss %>% mutate( Sex = factor(Sex), Year = factor(Year) )
# Update tibble
brfss <- brfss %>% mutate( Sex = factor(Sex), Year = factor(Year) )
brfss %>% filter(Sex == "Female", Year == "1990") %>% 
    select(Age, Weight, Height)
brfss %>%
    group_by(Sex, Year) %>% 
    summarize(
        AveAge = mean(Age, na.rm=TRUE),
        AveWeight = mean(Weight, na.rm=TRUE)
    )

# Import data
#-----------------------------------------------------#
pdata <- read_tsv("datasets/ALLphenoData.tsv")

# Statistical tests
#-----------------------------------------------------#
t.test(age ~ sex, pdata)
pdata %>% t.test(age ~ sex, data = .)
pdata %>% t.test(age ~ sex, data = ., var.equal = T)
# Develop function
t_test <- function(data, formula, ...) {
    t.test(formula, data, ...)
}
pdata %>% t_test(age ~ sex)

# Import data
#-----------------------------------------------------#
pdata_file <- file.choose() 
count_file <- file.choose()
pdata <- read_csv(pdata_file)
counts <- read_csv(count_file)

# Check data
#-----------------------------------------------------#
pdata <- 
    pdata %>% 
    select(Run, cell, dex)
pdata
eg <- counts[, 1:6] 
eg

# Combine data
#-----------------------------------------------------#
data <- left_join(pdata, eg)
data

# Transform data
#-----------------------------------------------------#
# Gather to transform the wide-format into long-format
tbl <- gather(data, "Gene", "Count", -(1:3))
tbl
tbl %>%
    group_by(Run) %>%
    summarize(lib_size = sum(Count))
tbl %>%
    group_by(Gene) %>%
    summarize(
        ave_count = mean(Count),
        ave_log_count = mean(log(1 + Count))
    )
counts_tbl <- gather(counts, "Gene", "Count", -Run)
data_tbl <- left_join(pdata, counts_tbl)
data_tbl
data_tbl %>%
    group_by(Run) %>%
    summarize(lib_size = sum(Count))
gene_summaries <-
    data_tbl %>%
    group_by(Gene) %>%
    summarize(
        ave_count = mean(Count),
        ave_log_count = mean(log(1 + Count))
    )
gene_summaries

# Visualize data
#-----------------------------------------------------#
ggplot(gene_summaries, aes(ave_log_count)) +
    geom_density()

# Session Info
#-----------------------------------------------------#
sessionInfo()