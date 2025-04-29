#### Consensus lineage detection functions ####

#' @importFrom igraph as.undirected cluster_louvain
#' @import proxy
#' @importFrom bezier bezier bezierArcLength
#' @importFrom cccd nng
#' @importFrom cluster pam
#' @importFrom stats heatmap
#'

# Organize scoring information and obtain membership matrices from the top-performing MSTs
get_top_scores <- function(all_scores, num_scores, percent) {
  s <- sapply(all_scores, "[")    # num_scores x k-possible matrix of data
  s <- unlist(sapply(s, "[")[2,]) # Extract scores
  p <- percent * num_scores * length(all_scores)
  top_scores <- order(s)[1:p]           # Identify  lowest scores
  clusn <- ifelse (top_scores %% num_scores == 0, top_scores / num_scores, floor(top_scores / num_scores) + 1)

  top_mm <- 1:p # top membership matrices
  for (i in 1:length(clusn)) {
    v <- all_scores[[clusn[i]]]     # number of scores / cluster
    index <- ifelse (top_scores[i] %% num_scores == 0, num_scores, top_scores[i] %% num_scores)
    top_mm[i] <- v[[index]][3]  # get membership matrix
  }

  top_data <- list('Top_MM' = top_mm, 'Top_Scores' = top_scores)
  return(top_data)
}


# Collects all membership matrices, grouped by Lineage 1, 2, etc.
get_lin_mem <- function(membership) {

  nl_detected <- unlist(lapply(membership, function(i) {
    x <- nrow(i)
    nl <- ifelse(is.null(x), 1, x)
    return(nl)
  })) # checks if cyclic/linear or multifurcating

  unl <- sort(unique(nl_detected))
  p2 <- max(nl_detected) # break 2

  b <- vector(mode='list', length=p2)[unl]
  b <- lapply(1:length(unl), function(i) {b[[i]] <- vector(mode='list', length=unl[i])})

  # Additional row of 0's to single-lineage membership "vector" --> matrix
  x <- lapply(1:length(membership), function(i)
    if (is.null(nrow(membership[[i]]))) {
      rbind(membership[[i]], rep(0, length(membership[[i]])))
    } else {
      membership[[i]]
    })

  b <- lapply(1:p2, function(i) {
    a <- x[which(nl_detected == i)] # Obtain lineages w/ "i" # of lineages
    for (j in 1:i) {
      df <- as.data.frame(sapply(a, "[", j,))
      colnames(df) <- seq_len(ncol(df))
      b[[j]] <- df
      b[[j]] <- t(as.data.frame(do.call(rbind, b[[j]])))
    }
    return(b)
  })

  all_membership <- unlist(b[unl], recursive = FALSE)
  z <- which(sapply(lapply(all_membership, nrow), is.null)) # remove unused lineages
  all_membership <- if(isempty(z)) {
    do.call(cbind, all_membership)
  } else { do.call(cbind, all_membership[-z]) }

  return(all_membership)
}


#' @title Determine consensus lineages in data.
#' @description Following the iterative MST inference across K clusters, this function
#' determines the optimal number of lineages based on the MST score and overlap
#' in paths that are detected by all MSTs.
#' @param retro_obj RETRO meta object
#' @param percent (optional) Fraction to select the top x% of high-performing MSTs based on score. Default is 0.05.
#' @param cutoff (optional) Membership cutoff for whether a cell can be assigned to a lineage based
#' on its mapping to that lineage across inferred MSTs. Default is 0.8 (cell must be mapped to lineage in ≥80% of selected MSTs)
#' @return This function returns \code{retro_obj} with updated information regarding the optimal MST for pseudotime inference, as well
#' as well as membership lists separating cells into their respective developmental lineages.
#' @export
#'
#'
# Determines # of Lineages (Heatmap) and Cell Membership to Each Lineage
get_num_lineages <- function(retro_obj, percent=0.05, cutoff=0.8) {

  all_k <- retro_obj@all_k # clustering information
  all_scores <- retro_obj@all_scores # scoring information
  threshold <- retro_obj@threshold # threshold used for terminal cluster specification
  num_scores = length(all_k[[1]]) # number of MSTs/clustering

  # Identify top scoring MSTs based on percent cutoff
  top_data <- get_top_scores(all_scores, num_scores, percent)
  membership <- top_data[[1]] # membership matrices
  top_scores = top_data[[2]]
  circular <-  isempty(unlist(lapply(membership, nrow))) # must all be ==1

  nl_detected <- unlist(lapply(membership, function(i) {
    x <- nrow(i)
    nl <- ifelse(is.null(x), 1, x)
    return(nl)
  })) # obtain lineages detected from top MSTs

  # (1) Cluster membership matrix
  if(!circular) {
    all_membership <- get_lin_mem(membership)

    # Jaccard distance on entire membership matrices for both lineages
    m <- as.matrix(proxy::dist(all_membership, by_rows = FALSE, method = "Jaccard"))
    rownames(m) <- colnames(m) <- 1:nrow(m)

    # Louvain clustering to identify # lineages
    knn <- igraph::as.undirected(nng(m, k=dim(m)[1]/max(nl_detected)))
    grouping <- cluster_louvain(knn)
    nl_knn <- length(grouping)

    heatmap(m)
    print(paste0('# of lineages (determined by Louvain clustering): ', nl_knn))

  } else {
    all_membership <- do.call(rbind, lapply(membership, function(y) if (is.vector(y)) t(y) else y))
    m = as.matrix(proxy::dist(all_membership, by_rows = FALSE, method = "Jaccard"))
    nl_knn <- nl <- 1
    print('All cells belong to 1 detected lineage.')
  }

  # (2) Cluster cells based on membership matrix
  factor = pam(m, nl_knn, diss=TRUE) # grouping for lineage membership matrices
  lin_membership <- lapply(1:nl_knn, function(i) {
    x <- t(all_membership[,factor$clustering==i])
    return(which(colMeans(x) >= cutoff))})

  # (3) Identify cells that map to lineage
  binary_mat = matrix(0, nrow=nrow(coordinates), ncol=length(lin_membership))
  for (j in seq_along(lin_membership)) {
    binary_mat[lin_membership[[j]], j] <- 1
  }

  # (4) Identify closest MST trajectory
  closest_mst <- sapply(1:length(membership), function(i) {
    membership_t <- t(membership[[i]])
    nl_t = ncol(membership_t)

    if(dim(membership_t)[1]==1) {
      membership_t = membership[[i]]
      nl_t = 1
    }
    compare_membership <- cbind(binary_mat, membership_t)
    mem_col = (nl_knn+1):ncol(compare_membership)

    compare_dist = as.matrix(proxy::dist(compare_membership, by_rows = FALSE, method = "Jaccard"))
    mn_dist = mean(compare_dist[1:nl_knn,mem_col]) + abs(nl_t-nl_knn)

    return(mn_dist)
  })

  # (5) Regenerate top-performing MST
  index = which.min(closest_mst) # based on membership matrix

  #### Secure specific clustering labels
  f <- top_scores[index]
  clusn.i <- ifelse (f %% num_scores == 0, f / num_scores, floor(f / num_scores) + 1) # which k
  iteration.i <- ifelse (f %% num_scores == 0, num_scores, f %% num_scores) # which iterations

  #### Obtain RETRO obj values to rerun MST calculation
  coordinates = retro_obj@coordinates
  kmedoids = all_k[[clusn.i]][[iteration.i]]

  retro_obj <- create_dMST(retro_obj, kmedoids)

  cell_MST = retro_obj@RETRO_MST
  id <- cell_MST$ID
  nl <- length(cell_MST$Lineages)
  centers <- lapply(1:nl, function(i) coordinates[id[unlist(cell_MST$Lineages[[i]])],])

  retro_obj@num_lineages <- nl # number of global lineages
  retro_obj@RETRO_MST <- cell_MST # optimal MST
  retro_obj@centroids <- centers # cell centroids used in optimal MST (for curve calculation)

  # Map cell MST to consensus lineage clustering
  # Calculate the average membership / lineage and keep cells above threshold
  m_matched_nl = t(do.call(rbind, membership[which(nl_detected==nl)]))
  m_matched_nl_dist = as.matrix(proxy::dist(m_matched_nl, by_rows = FALSE, method = "Jaccard"))

  factor_matched_nl = pam(m_matched_nl_dist, nl, diss=TRUE) # grouping for lineage membership matrices
  lin_membership <- lapply(1:nl, function(i) {
    x <- t(m_matched_nl[,factor_matched_nl$clustering==i])
    return(which(colMeans(x) >= cutoff))})

  cells_to_lin <- vector(mode='list', length=nl)
  cells_to_lin <- lapply(1:nl, function(i) {
    x <- t(m_matched_nl[,factor_matched_nl$clustering==i])
    return(coordinates[which(colMeans(x) >= cutoff),])})

  lineages <- cell_MST$Lineage
  if(nl > 1) {
    assign_lin <- lapply(lin_membership, function(l) {
      proportion <- lapply(lineages, function(lin) {
        prop <- length(which(l %in% id[lin]))
        return(prop)
      })
      return(proportion)
    })
    actual_order <- unlist(lapply(assign_lin, function(x) which.max(unlist(x))))
    actual_order <- order(actual_order)
    lin_membership <- lapply(actual_order, function(i) lin_membership[[i]])
    cells_to_lin <- lapply(actual_order, function(i) cells_to_lin[[i]])
  }

  retro_obj@cells_to_lin <- cells_to_lin
  retro_obj@lin_membership <- lin_membership

  return(retro_obj)
}



### Bézier spline fitting functions ####

# Extension for Bézier curve to cover all points
extend_centers <- function(centroids, p) {
  # Starting center extension
  v_start <- centroids[1,] - centroids[2,]
  root_center <- centroids[1,] + p * v_start

  # Ending center extension
  n <- nrow(centroids)
  v_terminal <- centroids[n,] - centroids[n-1,]
  terminal_center <- centroids[n,] + p * v_terminal

  extreme_centers <- rbind(root_center, terminal_center)

  return(extreme_centers)
}


#' @title Bézier spline fitting function
#' @description This function creates the Bézier spline curve to smooth out the
#' piecewise MST structure inferred by RETRO. It achieves this by stitching together cubic Bézier curves,
#' where the initial and final control points are from the nodes of the underlying MST (cluster "centroids").
#' The remaining two control points are computed internally and are not outputted by RETRO.
#' @param retro_obj List returned by RETRO \code{create_meta_obj()} containing information
#' about the optimal MST after scoring.
#' @param extension (optional) Numerical parameter for extending the curve beyond the
#' trajectory starting and terminal nodes. Default is 2.
#' @return List of curve information stored in \code{retro_obj@RETRO_Curve}.
#' The list elements are: (1) smoothed Bézier curves
#' <br>(2) MST centroids used to determine the curve
#' <br>(3) smoothed Bézier curves prior to extension
#' <br>(4) the vector of arc lengths between each set of centroids.
#' @export
#'
#'
# Contains concatenated Bézier curve functions, manages extension of Bézier curve
get_bezier_curve <- function(retro_obj, extension=2) {

  nl <- retro_obj@num_lineages

  bcurve_data <- lapply(1:nl, function(i) {
    centroids <- retro_obj@centroids[[i]]

    extreme_centers <- extend_centers(centroids, extension)
    centroids.2 <- rbind(extreme_centers[1,], centroids, extreme_centers[2,])

    ext_b <- stitch_bezier_curves(centroids.2)[[1]] # Obtains data points for EXTENDED Bezier curves
    path_b <- stitch_bezier_curves(centroids) # Obtain data points for regular Bezier curves
    arclength <- path_b[[3]]
    path_b <- path_b[[1]]

    sl <- norm(c(extreme_centers[1,1:2]-centroids[1,1:2]), type='2') # making this edit given myo5 compressed pt
    fl <- norm(c(extreme_centers[2,1:2]-centroids[nrow(centroids),1:2]), type='2')
    arclength <- c(sl, arclength, fl)

    return(list(ext_b, centroids, path_b, arclength))
  })

  retro_obj@RETRO_Curve <- bcurve_data
  return(retro_obj)

}

# Create Bézier curves between consecutive centroids inputted (in order)
stitch_bezier_curves <- function(S, inc=.001) {

  n <- nrow(S)
  d <- ncol(S)
  num_curves <- n-1

  # building 1-4 tridiagonal matrix
  M <- diag(4, num_curves)
  M[abs(row(M) - col(M)) == 1] <- 1
  M[1] <- 2
  M[num_curves * (num_curves - 1)] <- 2
  M[num_curves ** 2] <- 7

  # matrix of ctrl_points, P
  P_mat <- matrix(nrow=num_curves, ncol=d)
  for (i in 1:num_curves) {
    P_mat[i,] <- 2 * (2 * S[i,] + S[i+1,])
  }
  P_mat[1,] <- S[1,] + 2 * S[2,]
  P_mat[num_curves,] <- 8*S[n-1,] + S[n,]

  # solve for first unknown control point, A
  A <- solve(M, P_mat)

  # solve for second unknown control point, B
  B <- matrix(nrow=nrow(A), ncol=ncol(A))

  for (i in 1:(num_curves-1)) {
    B[i,] <- 2 * S[i+1,] - A[i+1,]
  }
  B[num_curves,] = (A[n-1,] + S[n,]) / 2

  ctrl_points <- lapply(1:num_curves, function(i) rbind(S[i,], A[i,], B[i,], S[i+1,]))

  # solve for cubic Bezier curves
  piecewise_bezier <- lapply(1:num_curves, function(i) { bezier(seq(0,1, inc), ctrl_points[[i]], deg=3) })
  arclength <- lapply(1:num_curves, function(i) { bezierArcLength(ctrl_points[[i]][,1:2], t1=0, t2=1, deg=3) })

  bcurve <- do.call(rbind, piecewise_bezier)
  arclength <- unlist(sapply(arclength, "[", 'arc.length'))

  return(list(bcurve, ctrl_points, arclength))
}

