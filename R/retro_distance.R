
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

# Weighted PC matrix
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

# Weight coordinates using realtime_dist functions (still load!!)
weight_coord <- function(scdata, weight) {
  # counts <- scdata$Counts
  time <- as.numeric(scdata$Time)
  pca <- scdata$PCA
  mat_list <- pc_distance_function(time, pca, weight)
  coordinates <- mat_list$Coordinates
  return(coordinates)
}

# Get average time point / cluster WITH TERMINAL CELL SUPERVISION
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
  threshold <- threshold * (k / max_k) # SCALE user-specified threshold
  
  ## TERMINAL CLUSTER SPECIFICATION
  # proportion of cells at final time-point (>x% to be terminal cluster)
  if(!isempty(terminal_cells)) { # based on labeling
    
    terminal_u <- lapply(seq(v), function(i) length(which(v[[i]]$id %in% unlist(terminal_cells))) / length(v[[i]]$id)) 
    terminal_id <- (names(v))[unlist(terminal_u) >= 0.10] # at least 10% labeled as terminal
    terminal_id <- as.numeric(terminal_id)
    
  } else { # based on time information
    terminal_t <- lapply(seq(v), function(i) length(which(v[[i]]$time==final_time)) / length(v[[i]]$time)) 
    terminal_id <- (names(v))[unlist(terminal_t) >= threshold] # at least 10% at final time-point
    terminal_id <- as.numeric(terminal_id)
  }
  
  # If thresholds are TOO RESTRICTIVE, any cluster w/ terminal cell = terminal cluster
  if(isempty(terminal_id)) {
    terminal_t <- lapply(seq(v), function(i) length(which(v[[i]]$time==final_time)) / length(time)) 
    terminal_id <- (names(v))[unlist(terminal_t) > 0] 
    terminal_id <- as.numeric(terminal_id)
  }
  
  ## STARTING CLUSTER SPECIFCATION
  if(!isempty(starting_cells)) { # based on labeling
    starting_u <- lapply(seq(v), function(i) length(which(v[[i]]$id %in% unlist(starting_cells))) / length(v[[i]]$id)) 
    start_id <- (names(v))[which.max(starting_u)] # at least 50% labeled as terminal
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

# Minimum Spanning Tree
create_dMST <- function(coordinates, kmedoids, time, 
                        terminal_cells=NULL, starting_cells=NULL,
                        threshold=0.10, max_k,
                        start=NULL, period=NULL) {
  
  # Identify Cluster Centroids
  clusterLabels <- kmedoids[["clustering"]]
  reclustered <- FALSE 
  
  if(!isempty(period)) {
    reclustering <- recluster(kmedoids, time, coordinates, period=period)
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
  
  # Set up for Minimum Spanning Tree
  g <- graph.adjacency(mat, mode = "undirected", weighted = TRUE)
  edge_list <- apply(get.edgelist(g), 2, as.numeric)
  weight <- E(g)$weight
  edge_list <- cbind(edge_list, weight)
  
  if(isempty(period)) { # Using timeStats/clusterLabels for non-cyclic data
    
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
      # n <- (1:length(nodes))[-node]
      node_b <- which(c(node.1, node.2) %in% b) # node B used to determine dt
      
      dt <- abs(av_time[node]-av_time[node_b]) # Difference in time between nodes
      return(dt)
    })
    dt <- ceiling(max(unlist(consec_dt))) / 2
  }
  
  ### PRUNING EDGE-WEIGHTED GRAPH ###
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
  
  # MINIMUM SPANNING TREE CALCULATION  
  g <- graph_from_edgelist(edge_list[,1:2], directed=FALSE)
  mst <- minimum.spanning.tree(g, weights=edge_list[,3])
  
  ### START/TERMINAL CLUSTER TROUBLESHOOTING
  if(!(start.clus %in% as_edgelist(mst))) {
    start.clus = as_edgelist(mst)[1]
  } 
  if (!prod(terminal.clus %in% as_edgelist(mst))) {
    terminal.clus = as.numeric(V(g))
  }
  paths <- all_simple_paths(mst, from=start.clus, to=terminal.clus)
  
  if(isempty(paths)) {
    start.clus <- as_edgelist(mst)[1]
    paths <- all_simple_paths(mst, from=start.clus)
  }
  ### 
  
  if(!isempty(filter_lineages(paths))) {paths <- filter_lineages(paths)} # redundant paths
  if(2 %in% lengths(paths) && length(paths)>1) { # v47
    paths = paths[-which(2 %in% lengths(paths))]
  }
  
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

