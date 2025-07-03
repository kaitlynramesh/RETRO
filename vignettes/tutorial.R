## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(RETRO)

## -----------------------------------------------------------------------------
library(Biobase) # data processing + extraction
library(ggplot2) # plotting
library(ggpubr) # plotting
library(doParallel)
library(parallel)

## -----------------------------------------------------------------------------
# Example 1: Bifurcating lineages
bifurcation_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-data/main/synthetic/scd_bifurcation.rds"
scd_bifurcation = readRDS(gzcon(url(bifurcation_url)))

# Time information (example 1 bifurcation)
time = phenoData(scd_bifurcation)@data[["time"]]
print(sort(unique(time)))

## -----------------------------------------------------------------------------
# PCA object
pca = experimentData(scd_bifurcation)@other[["PCA"]]

time_label = as.factor(time) # factor for coloring
ggplot(as.data.frame(pca$x), aes(x=PC1,y=PC2,colour=time_label)) + geom_point() + theme_bw()

