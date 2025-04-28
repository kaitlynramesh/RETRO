library(RETRO)

rm(list=ls())

load("~/KaitlynRRStudio/RETRO/data/scd_bifurcation.rda")

k_range <- seq(10,15)
lambda = 0
beta = 1
period <- NULL
time <- scd_bifurcation$Time

pca_obj = prcomp(scd_bifurcation$Counts, center=T)
scd_bifurcation$PCA = pca_obj

coordinates <- weight_coord(scd_bifurcation, weight=1)
coordinates_uw <- weight_coord(scd_bifurcation, weight=0)

retro_obj = set_RETRO_class(coordinates=coordinates,
                            coordinates_uw=coordinates_uw,
                            k_range = k_range,
                            time=time)

ncores <- parallel::detectCores()
registerDoParallel(cores=ncores)
retro_obj = scoring(retro_obj, k_range, num_scores = 10)

retro_obj <- get_num_lineages(retro_obj, percent=0.05)
retro_obj <- get_bezier_curve(retro_obj)
retro_obj <- get_mapped_cells(retro_obj)
retro_obj <- pseudotime_fit(retro_obj)

pseudotime <- retro_obj@pseudotime
curves = sapply(retro_obj@RETRO_Curve, "[", 3)

ggplot(as.data.frame(coordinates), aes(x=PC1, y=PC2, colour=pseudotime)) +
  geom_point() +
  scale_colour_continuous(low="blue", high="orange") +
  geom_point(data = as.data.frame(curves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(curves[[2]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()
