library(cluster)
library(pracma)
library(stats)
library(inflection)
library(parallel)
library(doParallel)
library(SimilarityMeasures)
library(Biobase)
library(cccd)
library(scales)

# Iteratively generates MSTs and scores them based on distance metrics
generate_paths <- function(kmedoids, time, coordinates, terminal_cells=NULL,
                           starting_cells=NULL, threshold=0.10,
                           max_k, start='Mode', period=NULL) {

  cell_MST <- create_dMST(coordinates, kmedoids, time, terminal_cells, starting_cells,
                          threshold, max_k, start, period)
  clusterLabels <- cell_MST$ClusterLabels

  points <- coordinates[,1:3] # Obtain 2D coordinates for PC projection
  cell_dist <- av_cell_dist(cell_MST=cell_MST, clusterLabels=clusterLabels, points=points, time=time) # Calculate distance score

  iteration <- list(cell_dist, kmedoids)

  return(iteration)
}


# Obtains average distance from cell to the lineages - generate_paths()
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
  cell_dist <- diag(matrix(cell_dist, nrow=nrow(edges),ncol=nrow(edges),byrow=TRUE))

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
    # index <- tidyr::replace_na(index, 0)
    index[is.na(index)] <- 0

    mem[i]<-ifelse(abs(index[1]-index[2])==1, 1,0)
  }
  return(mem)
}

### SCORING AND MST GENERATION ###

# Scoring uses foreach/dopar to parallelize WITH PENALIZATION
scoring <- function(coordinates, time, c1, k, num_scores, num_cells, terminal_cells=NULL, starting_cells=NULL, max_k, start='Average', threshold=0.10, period=period) {
  iteration <- foreach(i=1:num_scores, .combine='c') %dopar% {

    m <- kmnn_cluster(coordinates = c1, num_centers = k)
    iteration <- generate_paths(m, coordinates=coordinates, time=time,
                                terminal_cells=terminal_cells,
                                starting_cells=starting_cells, max_k=max_k,
                                start=start, threshold=threshold, period=period)
    list(iteration)
  }
  return(iteration)
}
