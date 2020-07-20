# Load tidyverse package
#-----------------------------------------------------#
library(readr)
library(tibble)
library(dplyr)
library(ggplot2)

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
t.test(age ~ sex, pdata)
pdata %>% t.test(age ~ sex, data = .)
pdata %>% t.test(age ~ sex, data = ., var.equal = T)
# Import data
#-----------------------------------------------------#


# Import data
#-----------------------------------------------------#