#### Distance and MST-related functions ####

# Determine number of PCs that make up 90% variance
num_pc <- function(pca) {
  n_pc <- 1:ncol(pca$x) # total pca

  eig <- pca$sdev ** 2
  contribution <- 1:length(n_pc)
  contribution[1] <- eig[1]
  for (i in 2:ncol(pca$x)) {
    contribution[i] <- (contribution[i-1] + eig[i])
  }
  contribution <- contribution / contribution[length(n_pc)]
  df <- data.frame(x=n_pc, y=contribution)

  curve <- lm(y~poly(x, 5), data=df) # degree 5 approx of contribution curve
  xval <-data.frame(x=seq(n_pc))
  approx <- predict(curve, newdata=xval)
  fitline <- data.frame(xval,approx)

  index <- which.min(abs(.9-fitline$approx)) # minimize diff between .9 and estimated values
  estimate <- fitline[index,] # find estimated coordinates where contribution reaches 90%

  e <- estimate$x
  if (e == 1) {e <- e+1}

  return (e)
}



#' @title Defining cell coordinates using gene expression and time
#' @description This function concatenates weighted time points to PCA coordinates.
#' @param time Vector of sampling time points corresponding to each cell.
#' @param pca List containing PCA information from \code{prcomp()}.
#' @param lambda A numerical value indicating contribution from PCA for scaling the \code{time} vector.
#' Default value is 1 (where PCA and time have an equal contribution). Increasing lambda will increase contribution from PCA.
#' @return List containing weighted coordinates, PCA coordinates weighted by their variance,
#' and the computed weight used to scale the time vector.
#'
pc_distance_function <- function(time,
                                 pca,
                                 lambda=1) {

  time = as.integer(factor(time))

  # Time distance matrix to calculate mean time between cells
  tmat <- dist(time) # Euclidean
  tmat <- as.matrix(tmat)
  tmat <- tmat**2 # difference squared
  mean_t <- mean(as.numeric(tmat))

  # PCs that contribute to 90% (default) of the variance in the data
  n_pc <- num_pc(pca)
  pc_sdev <- pca$sdev[1:n_pc]

  # PC values weighted by their variance
  pc_weighted <- apply(pca$x[,1:n_pc], 1, function(x) pc_sdev*x)

  pc_weighted <- t(pc_weighted)
  colnames(pc_weighted) <- colnames(pca$x[,1:n_pc])
  rownames(pc_weighted) <- rownames(pca$x[,1:n_pc])

  # PCA distance matrix to calculate median distance between PC values
  pc_mat <- as.matrix(dist(pc_weighted))
  pc_mat <- pc_mat**2
  median_pc <- median(as.numeric((pc_mat)))

  # PC factor that can equate the mean time distance to distance between PCs
  pc_factor <- median_pc / mean_t
  weight <- sqrt(pc_factor * lambda)

  # Weighted PC and time-coordinate system for cells
  coordinates <- cbind(pc_weighted, time * weight)

  return (list('Coordinates' = coordinates,
               'Weighted_PC' = pc_weighted,
               'Time_Weight' = weight))
}

#' @title Defining cell coordinates using gene expression and time
#' @description This function is a wrapper for the function \code{pc_distance_function()}.
#'
#' @param scdata A data frame containing experimental information and the PCA projection that will be used.
#' The vector of sampling time for each cell should be stored in scdata$Time and the
#' PCA object returned by \code{prcomp()} should be stored in scdata$PCA.
#' @param weight A numerical value indicating contribution from PCA for scaling the \code{time} vector.
#' Default value is 1 (where PCA and time have an equal contribution). Increasing lambda will increase contribution from PCA.
#' @return Returns the updated cell coordinates based on the PCA projection and sampling time.
#' @export
#'
# Weight coordinates using realtime_dist functions
weight_coord <- function(scdata, weight) {
  time <- as.numeric(scdata$Time)
  pca <- scdata$PCA
  mat_list <- pc_distance_function(time, pca, weight)
  coordinates <- mat_list$Coordinates
  return(coordinates)
}


#' @title Determing time information per cluster
#' @description This function infers the starting and terminal clusters based on time information.
#' Such information can also be user-specified if preferred.
#'
#' @param clusterLabels A membership vector for clustering the coordinates into groups.
#' @param time Vector of sampling time points corresponding to each cell.
#' @param terminal_cells (optional) List of cells belonging to the terminal ends of the trajectory.
#' Each list element can contain the cells corresponding to one terminal cell type.
#' @param starting_cells (optional) Vector of cells belonging to the starting end of the trajectory.
#' RETRO currently allows for only one starting node.
#' @param threshold (optional) Cutoff to classify whether clusters are terminal ends of the trajectory. Default is 0.1.
#' Note that this cutoff is scaled by number of clusters (K) / maximum number of clusters specified (Max_K).
#' @param max_k Maximum number of clusters from the user-specified range of K.
#' @return Returns consensus information about clusters such as their average and mode time (AvTime, Mode),
#' as well as the inferred initial and terminal clusters based on sampling time (Terminal_ID, Start_ID).
#' @export
#'
get_cluster_time <- function(clusterLabels, time, terminal_cells=NULL,
                             starting_cells=NULL, threshold=0.10, max_k) {

  id <- 1:length(clusterLabels)
  v <- as.data.frame(cbind(id, clusterLabels, time))
  v <- v[order(clusterLabels),]
  v <-split(v, v$clusterLabels)
  av_time <- 1:length(v)
  mode_time <- 1:length(v)
  final_time <- max(time) # assume that terminal clusters are present at LAST TIME

  k <- length(unique(clusterLabels))
  threshold <- threshold * (k / max_k) # SCALE threshold by maximum clusters

  ## TERMINAL CLUSTER SPECIFICATION
  # proportion of cells at final time-point (>x% to be terminal cluster)
  if(!isempty(terminal_cells)) { # based on labeling

    terminal_u <- lapply(seq(v), function(i) length(which(v[[i]]$id %in% unlist(terminal_cells))) / length(v[[i]]$id))
    terminal_id <- (names(v))[unlist(terminal_u) >= 0.10] # at least 10% labeled as terminal
    terminal_id <- as.numeric(terminal_id)

  } else { # based on time information
    terminal_t <- lapply(seq(v), function(i) length(which(v[[i]]$time==final_time)) / length(v[[i]]$time))
    terminal_id <- (names(v))[unlist(terminal_t) >= threshold] # at least x% at final time-point
    terminal_id <- as.numeric(terminal_id)
  }

  # If thresholds are TOO RESTRICTIVE, any cluster w/ terminal cell = terminal cluster
  if(isempty(terminal_id)) {
    terminal_t <- lapply(seq(v), function(i) length(which(v[[i]]$time==final_time)) / length(time))
    terminal_id <- (names(v))[unlist(terminal_t) > 0]
    terminal_id <- as.numeric(terminal_id)
  }

  ## STARTING CLUSTER SPECIFICATION
  if(!isempty(starting_cells)) { # based on labeling
    starting_u <- lapply(seq(v), function(i) length(which(v[[i]]$id %in% unlist(starting_cells))) / length(v[[i]]$id))
    start_id <- (names(v))[which.max(starting_u)]
    start_id <- as.numeric(start_id)
  } else {
    start_id <- NULL
  }

  for (i in 1:length(v)) {
    t <- v[[i]]$time
    o <- match(t, sort(unique(time))) # order
    av_time[i] <- mean(o) # pseudo-average time
    mode_time[i] <- pracma::Mode(o)
  }

  cluster_id <- as.double(names(v))

  return(list('AvTime'=av_time,
              'Mode'=mode_time,
              'Cluster_ID'=cluster_id,
              'Terminal_ID'=terminal_id,
              'Start_ID'=start_id))
}



#' @title Calculating the minimum spanning tree for trajectory inference
#' @description This function uses RETRO functions \code{recluster()} and \code{get_cluster_time()} to figure out
#' clustering and establish the initial and final nodes of the minimum spanning tree. It
#' depends on functions from the \code{igraph} package for most graph-based calculations.
#'
#' @param coordinates Matrix of cell coordinates with a time contribution
#' obtained using the \code{weight_coord()} function.
#' @param kmnn Cluster centroids and labels obtained from \code{kmnn_cluster()}.
#' @param time Vector of sampling time points corresponding to each cell.
#' @param start Character parameter that is either "Average" or "Mode." If "Average,"
#' the starting node is assigned to the cluster with the lowest average sampling time.
#' If "Mode," the starting node is assigned the cluster with the greatest representation of the lowest sampling time-point.
#' "Mode" is recommended for determining the starting node.
#' @param period Specifies the minimum difference in time where one can expect the cells
#' to return to an earlier state gene expression (make a cycle).
#' @return List of trajectory information, which is as follows: minimum spanning tree object (MST),
#' starting and terminal nodes (Start, Terminal), constituent lineages inferred (Lineages),
#' cluster labels for cell coordinates (clusterLabels), node IDs (ID),
#' MST edge cutoff based on time (Dt), and a logical value for whether the cells were reclustered (Reclus)
#' @export
#'
create_dMST <- function(coordinates, kmnn, time,
                        terminal_cells=NULL, starting_cells=NULL,
                        threshold=0.10, max_k,
                        start=NULL, period=NULL) {

  # Identify cluster centroids
  clusterLabels <- kmnn[["clustering"]]
  reclustered <- FALSE

  if(!isempty(period)) {
    reclustering <- recluster(kmnn, time, coordinates, period=period)
    clusterLabels <- reclustering[[1]]
    dt <- reclustering[[2]]
    period <- reclustering[[3]] # unchanged unless no separation
    reclustered <- ifelse(!isempty(period), TRUE, FALSE)
  }

  # Obtain average time point per cluster
  time_stats <- get_cluster_time(clusterLabels, time, terminal_cells,
                                 starting_cells, threshold=threshold, max_k)
  av_time <- time_stats[[1]]
  mode_time <- time_stats[[2]]

  id <- time_stats[[3]]
  medoids <- coordinates[id,]

  k <- length(id)
  nodes <- 1:k # node IDs

  # Identify starting cluster and terminal cluster
  if(!isempty(starting_cells)) {
    start.clus <- which(id %in% time_stats$Start_ID)
  } else {
    if(start=='Average') {
      start.clus <- which.min(av_time)
    } else if (start=='Mode') {
      start.clus <- which(mode_time == min(mode_time))
      start.clus <- ifelse(length(start.clus)>1, which.min(av_time), start.clus)
    }
  }
  terminal.clus <- which(id %in% time_stats$Terminal_ID)

  # Distance matrix for cluster centers
  mat <- as.matrix(dist(medoids))
  rownames(mat) <- colnames(mat) <- 1:length(id)

  # Set up for minimum spanning tree
  g <- graph.adjacency(mat, mode = "undirected", weighted = TRUE)
  edge_list <- apply(get.edgelist(g), 2, as.numeric)
  weight <- E(g)$weight
  edge_list <- cbind(edge_list, weight)

  if(isempty(period)) { # Determine "dt" cutoff specific to cyclic data

    # Determine difference in time between nodes to estimate dt
    consec_dt <- lapply(nodes, function(node) {
      node.1 <- which(edge_list[,2] %in% node) # not node IDs, row-numbers
      node.2 <- which(edge_list[,1] %in% node)

      dist_to_a <- edge_list[c(node.1, node.2), 3] # Distances incoming/outgoing from node A

      # Identify SECOND-CLOSEST node
      if(length(dist_to_a) == 1) {
        b <- c(node.1, node.2)[which.min(dist_to_a)] # if only 1 edge, *closest node used
      } else {
        b <- c(node.1, node.2)[order(dist_to_a)[2]]
      }
      node_b <- which(c(node.1, node.2) %in% b) # node B used to determine dt

      dt <- abs(av_time[node]-av_time[node_b]) # Difference in time between nodes
      return(dt)
    })
    dt <- ceiling(max(unlist(consec_dt))) / 2
  }

  # Determining MST from edge-weighted adjacency matrix
  # (1) Remove backward edges (∆t < 0) // should be functionally bi-directional
  pruned_edge_list <- lapply(nodes, function(node) {
    node_edges <- edge_list[c(which(edge_list[,1] %in% node),
                              which(edge_list[,2] %in% node)),]

    node_edges <- lapply(1:nrow(node_edges), function(i) {
      e <- as.numeric(node_edges[i,])
      dt_edge <- abs(av_time[e[1]]-av_time[e[2]]) # if ∆t is TOO far, remove edge
      e <- c(e, dt_edge)
      if(dt_edge >= dt) {
        e <- NULL
      }
      return(e)
    })
    node_edges <- do.call(rbind, node_edges)
  })
  pruned_edge_list <- do.call(rbind, pruned_edge_list)
  pruned_edge_list = unique(pruned_edge_list) # remove redundant edges
  if(!isempty(pruned_edge_list)) {
    edge_list = pruned_edge_list
  }

  # (2) Calculate Minimum Spanning Tree
  g <- graph_from_edgelist(edge_list[,1:2], directed=FALSE)
  mst <- minimum.spanning.tree(g, weights=edge_list[,3])

  # (3) Determine starting/terminal clusters (if incompatible w/ MST)
  if(!(start.clus %in% as_edgelist(mst))) {
    start.clus = as_edgelist(mst)[1]
  }
  if (!prod(terminal.clus %in% as_edgelist(mst))) {
    terminal.clus = as.numeric(V(g))
  }

  # (4) Identify paths from starting to terminal node in MST
  paths <- all_simple_paths(mst, from=start.clus, to=terminal.clus)
  if(isempty(paths)) {
    start.clus <- as_edgelist(mst)[1]
    paths <- all_simple_paths(mst, from=start.clus)
  }

  # (5) Remove redundant or 2-node paths
  if(!isempty(filter_lineages(paths))) {paths <- filter_lineages(paths)} # redundant paths
  if(2 %in% lengths(paths) && length(paths)>1) { # v47
    paths = paths[-which(2 %in% lengths(paths))]
  }

  # (6) Re-calculate terminal nodes from constituent paths of MST
  terminal.clus <- unlist(lapply(paths, function(x) tail(x,1)))

  path_list <- vector(mode='list', length=length(paths))
  names(path_list) <- paste(seq_len(length(paths)))
  for (i in 1:length(paths)) {
    path_list[[i]] <- paths[[i]]
  }

  return(list('MST' = mst,
              'Start' = id[start.clus],
              'Lineages' = path_list,
              'ID' = id,
              'Terminal' = id[terminal.clus],
              'ClusterLabels'=clusterLabels,
              'AvTime'=av_time,
              'Dt'=dt,
              'Reclus'=reclustered))
}


# Prevents repeated lineages from being outputted by MST
filter_lineages <- function(lineages) {

  nl <- length(lineages)
  l <- lapply(lineages, as.vector)
  x <- lapply(l, length)

  r <- lapply(seq(l), function(i) lapply(1:nl, function(j)
    sum((l[[j]] %in% l[[i]]))))
  r <- lapply(seq(r), function(i) which(r[[i]]==x[[i]]))
  repeated <- which(unlist(lapply(r, length)) > 1)
  lineages <- lineages[-repeated]

  return(lineages)
}

