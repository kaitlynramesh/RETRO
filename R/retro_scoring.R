#### MST generation and scoring functions ####

#' @importFrom parallel mclapply
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom pracma cross
#' @importFrom igraph get.edgelist
#'
#'
#' @title Infer and score MST
#' @description This function is a wrapper for key RETRO functions, namely
#' \code{create_dMST()} for MST inference and \code{av_cell_dist()} for MST scoring.
#' It is iteratively called in the function \code{scoring()} to generate MSTs for
#' coordinates clustered in a variety of ways.
#' @param retro_obj RETRO meta object.
#' @param kmedoids Cluster centroids and labels obtained from \code{kmnn_cluster()}.
#' @param max_k Maximum number of clusters from the user-specified range of K.
#' @return List containing the (1) MST scoring information from \code{av_cell_dist()}
#' and (2) clustering.
#' @export
#'
generate_paths <- function(retro_obj, kmedoids, max_k) {

  coordinates = retro_obj@coordinates
  time = retro_obj@time

  # determine cluster assignment and MST
  retro_obj = create_dMST(retro_obj, kmedoids)
  cell_MST = retro_obj@RETRO_MST
  clusterLabels = retro_obj@RETRO_MST[["ClusterLabels"]]

  points = coordinates[,1:3] # Obtain 2D coordinates for PC projection
  cell_dist = av_cell_dist(cell_MST=cell_MST, clusterLabels=clusterLabels, points=points, time=time) # Calculate distance score

  return(list(cell_dist, kmedoids))
}


#' @title Computes a distance-based score for each MST
#' @description This function projects each cell to each lineage of the MST to gauge how well the
#' MST covers the low-dimensional projection of data. It includes penalties based on whether all of the
#' nodes are used, and if the initial and terminal nodes are consistent with the experimental sampling time.
#'
#' @param cell_MST List returned by \code{create_dMST()} containing MST-specific information, such as
#' MST structure, starting/terminal nodes, etc.
#' @param clusterLabels A membership vector for clustering the coordinates into groups.
#' @param points Coordinates for each cell. For low-dimensional projections from PCA, these coordinates use
#' the first 3 PC's and will have no time contribution. Note that if UMAP is used, there will be a time contribution as the first
#' three dimensions are UMAP1, UMAP2, and the weighted time vector.
#' @param time Vector of sampling time points corresponding to each cell.
#' @return List containing the (1) number of lineages, (2) score, and a (3) binary membership matrix indicating
#' which lineage the cell is closest to. This membership matrix will be used for detecting global lineages based on
#' all MSTs.
#'
av_cell_dist <- function(cell_MST, clusterLabels, points, time) {

  mst <- cell_MST$MST
  id <- cell_MST$ID
  lineages <- cell_MST$Lineages
  av_time <- cell_MST$AvTime
  reclustered <- cell_MST$Reclus # to check overlapping paths
  medoids <- points[id,]

  edges <- as.matrix(igraph::get.edgelist(mst)) # MST Edge Connection

  # Number of lineages
  nl <- length(lineages)

  # Distance between cell and its closest edge
  cell_dist <- apply(points, get_distance, MARGIN=1, edges=edges, center_coord=medoids)
  dist_score <- mean(cell_dist[1,]) # Average minimum distance

  # Identifies closest edge to cell
  closest_edge <- cell_dist[2,]
  closest_edge <- edges[closest_edge,]

  # Cell is assigned to lineage of closest edge
  cell_membership <- apply(closest_edge, is.member, MARGIN=1,lineages=lineages)

  # Utilization Penalty
  u <- length(unique(unlist(lineages))) # number of utilized nodess

  # Starting and terminal node penalties
  term_time <- av_time[id %in% cell_MST$Terminal] # time @ terminal node
  start_time <- av_time[id %in% cell_MST$Start] # time @ starting node
  tp <- max(abs(term_time - length(unique(time)))) # terminal penalty
  sp <-  max(abs(start_time - 1)) # starting penalty

  # Calculate penalized score
  penalized_score <- dist_score + (nrow(medoids) - u) + (tp+sp)

  return(list(nl, penalized_score, cell_membership))
}

# Determines minimum distance of cell to lineages and the corresponding edge
get_distance <- function(point, edges, center_coord) {

  dimensions = c(1,2,3)
  n1 = dimensions[1]
  n2 = dimensions[2]
  n3 = dimensions[3]

  x <- as.numeric(point[n1])
  y <- as.numeric(point[n2])
  z <- as.numeric(point[n3])

  cell_dist <- 1:nrow(edges) # Initialize vectors

  check_alpha <- function(i) {

    c1 <- as.numeric(edges[i,1]) # Identify start/end clusters per segment
    c2 <- as.numeric(edges[i,2])

    # Obtain coordinates for the segment endpoints
    x1 <- as.numeric(center_coord[c1,])[1]
    y1 <- as.numeric(center_coord[c1,])[2]
    x2 <- as.numeric(center_coord[c2,])[1]
    y2 <- as.numeric(center_coord[c2,])[2]
    p1 <- as.numeric(center_coord[c1,])
    p2 <- as.numeric(center_coord[c2,])

    # Calculate alpha
    a <- -((y-y2)*(y2-y1) + (x-x2)*(x2-x1)) / ((x1-x2)**2 + (y1-y2)**2)

    # Calculating 3D distance between points
    v1 <- center_coord[c1,dimensions]-point
    v2 <- center_coord[c2,dimensions]-point
    d1 <- norm(v1, type='2') # Distance 1 (a < 0)
    d2 <- norm(v2, type='2') # Distance 2 (a > 1)

    n <- p1[dimensions] - p2[dimensions]
    min_dist <- norm(pracma::cross(v1, n), type='2') / norm(n, type='2') # Minimum distance (perpendicular)

    if (a < 0 || a > 1) { # Assign distance depending on alpha-value
      cell_dist[i] <- min(d1, d2) # Distance to closest endpoint if outside of segment
    } else {
      cell_dist[i] <- min_dist # Minimum distance if within segment
    }
    return(cell_dist)
  }

  i <- 1:nrow(edges)
  cell_dist <- unlist((lapply(i, check_alpha))) # Find distance between point and all edges
  cell_dist <- diag(matrix(cell_dist, nrow=nrow(edges),ncol=nrow(edges), byrow=TRUE))

  closest_edge <- match(min(cell_dist), cell_dist)

  cell_dist <- min(cell_dist) # Return smallest distance between point and edge
  cell_membership <- c(cell_dist, closest_edge)

  return(cell_membership)
}

# Returns 1,0 for membership of cell to a lineage
is.member <- function(edge, lineages) {

  edge <- as.numeric(edge)
  n <- length(lineages)
  mem <- 1:n
  for (i in 1:n) {
    lineage <- unlist(lineages[[i]], use.names = FALSE)
    index <- match(edge, lineage)
    index[is.na(index)] <- 0

    mem[i]<-ifelse(abs(index[1]-index[2])==1, 1,0)
  }
  return(mem)
}


#' @title Parallelization of MST generation and scoring
#' @description This function is a wrapper for the \code{generate_paths()} function
#' to parallelize iterative MST generation and scoring across the range of specified K clusters.
#' It is recommended to run \code{scoring()} in an lapply() to contain all
#' iterative clustering and scoring data in one .rda file.
#'
#' @param retro_obj RETRO meta object.
#' @param num_scores Number of MSTs (scored MSTs) per clustering. Default is 100.
#' @param k_range Numeric vector for the range of K clusters that will be used to partition the coordinates
#' @return Lists in \code{retro_obj@all_scores} containing MST scoring information from \code{av_cell_dist()}
#' and in \code{retro_obj@k} for storing results from iterative clustering.
#' @export
#'
scoring <- function(retro_obj, k_range, num_scores=100) {

  ncores <- parallel::detectCores()
  registerDoParallel(cores=ncores)

  results <- lapply(k_range, function(k) {
    iterations <- foreach(i=1:num_scores, .combine='c') %dopar% {

      coordinates = retro_obj@coordinates
      coordinates_uw = retro_obj@coordinates_uw
      time = retro_obj@time
      terminal_cells = retro_obj@terminal_cells
      starting_cells = retro_obj@starting_cells
      start=retro_obj@start
      threshold=retro_obj@threshold
      period=retro_obj@period

      m <- kmnn_cluster(coordinates = coordinates_uw, num_centers = k)
      iteration <- generate_paths(retro_obj, kmedoids = m)
      list(iteration)
    }
  })

  # Extract scoring results
  all_scores <- lapply(seq_along(k_range), function(i) {
    l <- lapply(1:num_scores, function(j) results[[i]][[j]][[1]])
    return(l)})

  # Extract clustering results
  all_k <- lapply(seq_along(k_range), function(i) {
    l <- lapply(1:num_scores, function(j) results[[i]][[j]][[2]])
    return(l)})

  retro_obj@all_k = all_k
  retro_obj@all_scores = all_scores

  return(retro_obj)
}

