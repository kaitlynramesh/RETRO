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
library(ggpubr) # plotting
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

## -----------------------------------------------------------------------------
if (!interactive()) {
  foreach::registerDoSEQ()
} else {
  ncores = parallel::detectCores()
  registerDoParallel(cores = ncores)
}
retro_obj = scoring(retro_obj, k_range, num_scores)

boxplot_scoring(retro_obj) # visualization distribution of score over K clusters

## -----------------------------------------------------------------------------
retro_obj <- get_num_lineages(retro_obj, percent=0.1) 

## -----------------------------------------------------------------------------
retro_obj <- get_bezier_curve(retro_obj)
retro_obj <- get_mapped_cells(retro_obj) # map any left over cells to lineage

## -----------------------------------------------------------------------------
retro_obj = pseudotime_fit(retro_obj)

cyclic_lineage = retro_obj@pseudotime_graph[[1]]
cyclic_lineage

## -----------------------------------------------------------------------------
coordinates = retro_obj@coordinates
pseudotime = retro_obj@pseudotime
curve = retro_obj@RETRO_Curve[[1]][["Full_Curve"]] # curve data points

ggplot(as.data.frame(coordinates), aes(x=PC1, y=PC2, colour=pseudotime)) +
  geom_point() +
  scale_colour_continuous(low="blue", high="orange") +
  geom_point(data = as.data.frame(curve[,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()

## -----------------------------------------------------------------------------
pseudotime_density(time, pseudotime)

## -----------------------------------------------------------------------------
sessionInfo()

