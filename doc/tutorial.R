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

## -----------------------------------------------------------------------------
lambda = 1

coordinates <- weight_coord(scd_bifurcation, weight = lambda)
coordinates_uw <- weight_coord(scd_bifurcation, weight = 0)

## -----------------------------------------------------------------------------
k_range = seq(9,14) # set range of clusters 
num_scores = 20

retro_obj = set_RETRO_class(coordinates=coordinates,
                            coordinates_uw=coordinates_uw,
                            k_range = k_range,
                            time=time)

if (!interactive()) {
  foreach::registerDoSEQ()
} else {
  ncores = parallel::detectCores()
  registerDoParallel(cores = ncores)
}
retro_obj = scoring(retro_obj, k_range, num_scores)

## -----------------------------------------------------------------------------
boxplot_scoring(retro_obj)

## -----------------------------------------------------------------------------
retro_obj <- get_num_lineages(retro_obj, percent=0.05) 

## -----------------------------------------------------------------------------
retro_obj <- get_bezier_curve(retro_obj)

## -----------------------------------------------------------------------------
retro_obj <- get_mapped_cells(retro_obj)

cells_to_lin1 = retro_obj@cells_to_lin[[1]]
cells_to_lin2 = retro_obj@cells_to_lin[[2]]

lineage_1 = ggplot(as.data.frame(coordinates), aes(x=PC1, y=PC2)) +
  geom_point() +
  geom_point(data = as.data.frame(cells_to_lin1[,1:2]), aes(x = PC1, y = PC2), colour='#344AB0') +
  theme_bw()

lineage_2 = ggplot(as.data.frame(coordinates), aes(x=PC1, y=PC2)) +
  geom_point() +
  geom_point(data = as.data.frame(cells_to_lin2[,1:2]), aes(x = PC1, y = PC2), colour='#FF69B4') +
  theme_bw()

ggarrange(lineage_1, lineage_2, ncol=2)

## -----------------------------------------------------------------------------
retro_obj = pseudotime_fit(retro_obj)

lin_1 = retro_obj@pseudotime_graph[[1]]
lin_2 = retro_obj@pseudotime_graph[[2]]

ggarrange(lin_1, lin_2, ncol=2)

## -----------------------------------------------------------------------------
coordinates = retro_obj@coordinates
pseudotime = retro_obj@pseudotime
curves = sapply(retro_obj@RETRO_Curve, "[", 1)

ggplot(as.data.frame(coordinates), aes(x=PC1, y=PC2, colour=pseudotime)) +
  geom_point() +
  scale_colour_continuous(low="blue", high="orange") +
  geom_point(data = as.data.frame(curves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(curves[[2]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()

## -----------------------------------------------------------------------------
pseudotime_density(time, pseudotime)

## -----------------------------------------------------------------------------
sessionInfo()

