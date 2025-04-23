#### Visualization and plotting functions ####

# Creates data frame of center coordinates (df) in order of MST
get_mst_coord <- function(lineage, center_coord) {
  mstc <- data.frame()               # Start w/ empty data frame
  for (i in 1:length(lineage)) {
    x1 <- center_coord[lineage[[i]], 1] # Finds (x,y) for each center
    y1 <- center_coord[lineage[[i]], 2]
    c <- c(x1, y1)
    mstc <- rbind(mstc, c)
  }
  colnames(mstc) <- c('x', 'y')
  return(mstc) # Each lineage now is a data frame of cluster center coordinates
}

# Creates line segments to connect centers for plotting MST
connect_clus <- function(b, mstc) {
  for (i in 1:length(mstc)) { # Operates on each lineage identified
    for (j in 1:(length(mstc[[i]][["x"]])-1)) { # excessive indexing due to "list" object
      b <- b + geom_segment(x=mstc[[i]][["x"]][[j]],
                            y=mstc[[i]][["y"]][[j]],
                            xend=mstc[[i]][["x"]][[j+1]],
                            yend=mstc[[i]][["y"]][[j+1]])
    }
  }
  return(b) # Graph of MST
}

#' @import ggplot2
#' @import dplyr
#' @import tibble
#'
#' @title Boxplots to visualize score distribution
#' @description This function plots a boxplot per K number of clusters to look at
#' the distribution of scores of MSTs inferred when coordinates are partitioned into K clusters.
#' This plot in part explains the choice of K used in the optimal MST for curve fitting.
#'
#' @param all_scores List of MSTs obtained from the function \code{scoring()} that
#' is organized by the corresponding K clusters used. The length of the list should be
#' number of K clusters used, \code{length(k)}, and each list element should be the length \code{num_scores}.
#' The list contains the scores and structures of each MST determined.
#' @param num_scores Number of MSTs (scored MSTs) per clustering. Default is 100.
#' @param k Range of clusters used for RETRO pseudotime analysis.
#' @return A ggplot2 object for a set of boxplots of scores per K clusters
#' to show the effect of clustering on MST inference.
#' @export
#'
boxplot_scoring <- function(all_scores, num_scores, k) {

  kscoring <- lapply(all_scores, function(x) sapply(x, "[[", 2))
  score_df <- unlist(rbind(as.data.frame(kscoring)))
  num_k <- unlist(lapply(1:length(k), function(i) rep(k[i], num_scores)))
  score_df <- data.frame(cbind(unlist(score_df), num_k))
  rownames(score_df) <- 1:(num_scores * length(k))
  colnames(score_df) <- c('Score', 'K')
  score_df$K <- as.factor(score_df$K)
  ggplot(score_df, aes(x=K, y=Score)) +
    geom_boxplot(fill='steelblue1') + ggtitle('MST Structure Error from Iterative K-Medoids Clustering') +
    theme_bw() +
    theme(axis.text=element_text(size=10),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size=20))

}


# plot MST
test_MST <- function(paths, medoids, coordinates, clusterLabels) {

  center_coord <- data.frame(x=as.vector(medoids[,1]), y=as.vector(medoids[,2]))
  colnames(center_coord) = c("x", "y")
  colnames(coordinates) <- paste0('V', 1:ncol(coordinates))
  title <- paste0('Num Paths: ', length(paths))

  b <- ggplot(as.data.frame(center_coord), aes(x = x, y = y)) +
    geom_point(data = as.data.frame(coordinates)[,1:2], aes(x = V1, y = V2, color = as.factor(clusterLabels)), size = 2) +
    geom_point() +
    ggtitle(label=title) +
    theme_classic()

  lineages <- vector(mode='list', length=length(paths))
  names(lineages) <- paste(seq_len(length(paths)))
  for (i in 1:length(paths)) {
    lineages[[i]] <- paths[[i]]
  }

  lineages <- lapply(lineages, as.double)
  mstc <- lapply(lineages, get_mst_coord, center_coord=center_coord) # Identifies MST coordinates per lineage
  connect_clus(b, mstc)
}


#' @importFrom tidyr nest unnest
#' @importFrom purrr map2
#' @importFrom RColorBrewer brewer.pal

#' @title Pseudotime density plot
#' @description This function creates a density plot of pseudotime that is grouped
#' by discrete experimental time labels. This plot gives insight into whether the inferred
#' pseudotime is consistent with the progression of real time.
#' @param time Vector of experimental sampling time assigned to each cell.
#' @param pseudotime Vector of pseudotime estimates assigned to each cell.
#' @param bw (optional) Density plot bandwidth. Default at 0.5.
#' @return A ggplot2 object showing the normalized density of pseudotime values grouped by sampling time.
#' @export
#'
psupertime_density <- function(time, pseudotime, bw=0.5) {

  x <- seq(length(time)) # number of cells
  time <- as.numeric(time)
  rownames(x) <- paste0('cell-', 1:nrow(x)) # cell ids
  label_input = factor(time) # sampling time label

  proj_dt <- list(cell_id = rownames(x),
                  psuper = pseudotime,
                  label_input = label_input)
  proj_dt = as.data.frame(proj_dt)
  n_labels <- length(unique(label_input))

  if (n_labels <= 11) {
    col_vals 	= rev(brewer.pal(n_labels, name='RdBu'))
  } else {
    col_pal <- rev(brewer.pal(11, 'RdBu'))
    col_vals <- rev(colorRampPalette(col_pal)(n_labels))
  }

  bw_list = proj_dt %>%
    group_by(label_input) %>%
    summarize(bw = density(psuper)$bw)
  bw_list$bw = rep(bw, length(unique(label_input)))

  # (1) Table of custom bandwidths
  bw_list <- proj_dt %>%
    group_by(label_input) %>%
    summarize(bw = bw, .groups = "drop")

  # (2) Compute density per group from bandwidth list
  dens_df <- proj_dt %>%
    group_by(label_input) %>%
    nest() %>%
    left_join(bw_list, by = "label_input") %>%
    mutate(
      density_data = map2(data, bw, ~ {
        d <- density(.x$psuper, bw = .y)
        tibble(x = d$x, y = d$y / mean(d$y), label_input = .x$label_input[1])
      })
    ) %>%
    select(label_input, density_data) %>%
    unnest(density_data)

  bw_list_new = as.data.frame(dens_df) %>%
    group_by(label_input) %>%
    summarize(bw = density(x)$bw)
  # print(bw_list_new)

  ggplot(dens_df, aes(x = x, y = y, color = label_input, fill = label_input)) +
    geom_line() +
    geom_area(alpha = 0.4, position = "identity") +
    scale_colour_manual(values = col_vals) +
    scale_fill_manual(values = col_vals) +
    labs(
      x = "Pseudotime",
      y = "Density",
      color = "Time",
      fill = "Time"
    )

}
