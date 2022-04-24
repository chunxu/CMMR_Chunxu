# This document illustrates the main features of the DirichletMultinomial package, and in the process replicates key tables and gures from [1].
# [1] I. Holmes, K. Harris, and C. Quince. Dirichlet multinomial mixtures: Generative models for microbial metagenomics. PLoS ONE, 7(2):e30126, 02 2012.
# 

#install package DirichletMultinomial

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DirichletMultinomial")

# #install package xtable
# install.packages("xtable", repos="http://R-Forge.R-project.org")

#library
library(DirichletMultinomial)
library(lattice)
library(xtable)
library(parallel)
library("readxl")
library(rlang)
#setup format, etc.
options(width=70, digits=2)
full <- FALSE
.qualitative <- DirichletMultinomial:::.qualitative
dev.off <- function(...) invisible(grDevices::dev.off(...))

#load data
fl <- system.file(package="DirichletMultinomial", "extdata",
                  "MILES_V1.csv")
count <- t(as.matrix(read.csv(fl, row.names=1)))
count[1:5, 1:3]

#taxon-counts: distribution of reads from each taxon, on a log scale.
cnts <- log10(colSums(count))
getwd()
pdf("MILES_V1_taxon-counts.pdf")
densityplot(cnts, xlim=range(cnts),
            xlab="Taxon representation (log 10 count)")
dev.off()

#fit model
# if (full) {
#   fit <- mclapply(1:7, dmn, count=count, verbose=TRUE)
#   save(fit, file=file.path(tempdir(), "fit.rda"))
# } else data(fit)
# fit[[4]]
fit <- mclapply(1:7, dmn, count=count, verbose=TRUE)
fit[[4]]
fit[[3]] #minimum

# min-laplace
lplc <- sapply(fit, laplace)
pdf("MILES_V1_min-laplace.pdf")
plot(lplc, type="b", xlab="Number of Dirichlet Components",
     ylab="Model Fit")
dev.off()
(best <- fit[[which.min(lplc)]])

#mix-weight
mixturewt(best)
head(mixture(best), 3)

#fitted
pdf("MILES_V1_fitted.pdf")
splom(log(fitted(best)))
dev.off()

#posterior-mean-diff
p0 <- fitted(fit[[1]], scale=TRUE)     # scale by theta
p3 <- fitted(best, scale=TRUE)
colnames(p3) <- paste("m", 1:3, sep="")
(meandiff <- colSums(abs(p3 - as.vector(p0))))
sum(meandiff)

#table-1
diff <- rowSums(abs(p3 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 10)

#xtable
xtbl <- xtable(df,
               caption="Taxonomic contributions (10 largest) to Dirichlet components.",
               label="tab:meandiff", align="lcccccc")
print(xtbl, hline.after=0, caption.placement="top")

#heatmap-similarity
pdf("MILES_V1_heatmap1.pdf")
heatmapdmn(count, fit[[1]], best, 30)
dev.off()









