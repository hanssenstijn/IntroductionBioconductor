# Select file
#-----------------------------------------------------#
fl <- file.choose()

# Read data
#-----------------------------------------------------#
symgo <- read.csv(fl, row.names=1, stringsAsFactors=FALSE)

# Inspect data
#-----------------------------------------------------#
head(symgo)
dim(symgo)
length(unique(symgo$SYMBOL))

# Split data
#-----------------------------------------------------#
# split = divides data in vector x into the groups defined by f
go2sym <- split(symgo$SYMBOL, symgo$GO)

# Check length data
#-----------------------------------------------------#
len1 <- sapply(go2sym, length)  
len2 <- lengths(go2sym)
identical(len1, len2)
# Splits data into subsets, computes summary/statistics
len3 <- aggregate(SYMBOL ~ GO, symgo, length)
head(len3)
# Other way around
head(aggregate(GO ~ SYMBOL, symgo, length))
head(aggregate(SYMBOL ~ GO, symgo, c))

# Functions
#-----------------------------------------------------#
uidfun  <- function(x) {
    unique(tolower(x))
}
head(aggregate(SYMBOL ~ GO , symgo, uidfun))
# 'Anonymous' function
head(aggregate(SYMBOL ~ GO, symgo, function(x) {
    unique(tolower(x))
}))

# Phenotype data
#-----------------------------------------------------#
fname <- file.choose()
stopifnot(file.exists(fname))
pdata <- read.delim(fname)
class(pdata)

# Inspect data
#-----------------------------------------------------#
colnames(pdata)
dim(pdata)
head(pdata)
table(pdata$sex)
summary(pdata$cyto.normal)
pdata[1:5, 3:4]
pdata[1:5, ]
head(pdata[, 3:5])
tail(pdata[, 3:5], 3)
head(pdata$age)
head(pdata$sex)
head(pdata[pdata$age > 21,])
idx <- pdata$sex == "F" & pdata$age > 40
table(idx)
dim(pdata[idx,])

# Subset data
#-----------------------------------------------------#
bcrabl <- pdata[pdata$mol.biol %in% c("BCR/ABL", "NEG"),]
# Drop old factors, not used factors anymore
bcrabl$mol.biol <- factor(bcrabl$mol.biol)
table(bcrabl$BT)

# Analyze data
#-----------------------------------------------------#
aggregate(age ~ mol.biol + sex, bcrabl, mean)
t.test(age ~ mol.biol, bcrabl)

# Plot data
#-----------------------------------------------------#
boxplot(age ~ mol.biol, bcrabl)

# Weighty data
#-----------------------------------------------------#
wname <- file.choose()
stopifnot(file.exists(wname))
brfss <- read.csv(wname)

# Visualize data
#-----------------------------------------------------#
plot(sqrt(Weight) ~ Height, brfss, main="All Years, Both Sexes")

# Subset data
#-----------------------------------------------------#
brfss2010 <- brfss[brfss$Year == "2010", ]

# Visualize data
#-----------------------------------------------------#
opar <- par(mfcol=c(1, 2))
plot(sqrt(Weight) ~ Height, brfss2010[brfss2010$Sex == "Female", ],
     main="2010, Female")
plot(sqrt(Weight) ~ Height, brfss2010[brfss2010$Sex == "Male", ],
     main="2010, Male")
par(opar)     
# Load package
library(ggplot2)
qplot(Height, sqrt(Weight), data=brfss)
ggplot(brfss, aes(x=Height, y=sqrt(Weight))) +
    geom_point()

# Checking missing data
#-----------------------------------------------------#
sum(is.na(brfss$Height))
sum(is.na(brfss$Weight))
# Remove missing data
drop <- is.na(brfss$Height) | is.na(brfss$Weight)
sum(drop)
brfss <- brfss[!drop,]

# Re-plot data
#-----------------------------------------------------#
qplot(Height, sqrt(Weight), data=brfss) +
    ylab("Square root of Weight") + 
    ggtitle("All Years, Both Sexes")
ggplot(brfss, aes(x=Height, y=sqrt(Weight), color=Sex)) + 
    geom_point()
ggplot(brfss, aes(x=Height, y = sqrt(Weight), color=Sex, shape=Sex)) + 
    geom_point()
ggplot(brfss, aes(x=Height, y = sqrt(Weight), color=Sex)) + 
    geom_point() +
    facet_grid(Sex ~ .)
# Create subset
brfss2010 <- brfss[brfss$Year == "2010", ]
ggplot(brfss2010, aes(x=sqrt(Weight), fill=Sex)) +
    geom_density(alpha=.25)