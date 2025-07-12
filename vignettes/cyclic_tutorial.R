## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(RETRO)

# Additional packages for vignette 
library(Biobase) # data storage + processing
library(ggplot2) # plotting
library(doParallel)
library(parallel)

## -----------------------------------------------------------------------------
cyclic_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-data/main/synthetic/scd_twocycles.rds"
scd_twocycles = readRDS(gzcon(url(cyclic_url)))

## -----------------------------------------------------------------------------
# PCA analysis
pca = experimentData(scd_twocycles)@other$PCA

time = phenoData(scd_twocycles)@data[["time"]] 
time_label = as.factor(time) # factor for coloring
ggplot(as.data.frame(pca$x), aes(x=PC1,y=PC2,colour=time_label)) + geom_point() + theme_bw()

## -----------------------------------------------------------------------------
lambda = 1

coordinates <- weight_coord(scd_twocycles, weight = lambda)
coordinates_uw <- weight_coord(scd_twocycles, weight = 0)

## -----------------------------------------------------------------------------
k_range = seq(9,14) # set range of clusters 
num_scores = 20
period = 4

retro_obj = set_RETRO_class(coordinates=coordinates,
                            coordinates_uw=coordinates_uw,
                            k_range=k_range,
                            time=time,
                            period=period)

