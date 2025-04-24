#### Clustering functions ####

#' @importFrom stats dist kmeans density
#'
#'

find_nearest_cell <- function(centroid, coordinates) {
  dims = c(1,2,ncol(coordinates))
  centroid = centroid[dims]

  distances = apply(coordinates[,dims], 1, function(row) {
    return(as.numeric(stats::dist(rbind(row, centroid)))) # use all dimensions!!
  })
  nearest_idx <- which.min(distances)
  return(nearest_idx)
}

#' @title Clustering strategy adapted from K-means.
#' @description This function is an adaptation of k-means that maps K-means centroids
#' in each cluster to the nearest cell from the original data.
#'
#' @param num_centers Determines number of clusters.
#' @param coordinates Matrix of cell coordinates with a time contribution
#' obtained using the \code{weight_coord()} function.
#' @return A list containing the new cluster centroids \code{centers} and
#' a vector of labels \code{clusterLabels} grouping the coordinates.
#' @export
#'
#'
kmnn_cluster <- function(coordinates, num_centers) {
  res = kmeans(coordinates, centers = num_centers)
  centroids = res$centers
  clusterLabels = as.factor(res$cluster)

  nearest_cells <- sapply(1:num_centers, function(i) find_nearest_cell(centroids[i,], coordinates))
  nearest_cells = as.numeric(nearest_cells)

  levels(clusterLabels) = nearest_cells
  clusterLabels = as.numeric(as.vector(clusterLabels))

  return(list("centers"=nearest_cells,
              "clustering"=clusterLabels))
}


#' @title Re-cluster cells given a period
#' @description This function re-clusters the cells if there is periodicity in their gene expression,
#' causing clusters to contain multiple time points.
#'
#' @param kmnn Initial cluster centroids and labels from \code{kmnn_cluster()}.
#' @param period Specifies the minimum difference in time where one can expect the cells
#' to return to an earlier state gene expression (make a cycle).
#' @return A list containing an updated membership vector \code{clusterLabels} following re-clustering,
#' the average difference in time between cycles, \code{cycle_dt}, and the \code{period}.
#' @export
#'
recluster <- function(kmnn, time, coordinates, period) {

  clusterLabels <- kmnn[["clustering"]]
  bandwidth <- density(time)$bw
  id <- 1:length(clusterLabels)
  v <- as.data.frame(cbind(id, clusterLabels, time))
  v <- v[order(clusterLabels),]
  v <-split(v, v$clusterLabels)

  cuts <- vector(mode='list', length=length(v))
  prev_peaks <- vector(mode='list', length=length(v))
  peak_storage <- vector(mode='list', length=length(v))

  new_peaks <- prev_peaks
  names(prev_peaks) <- 1:length(v)

  for (i in 1:length(v)) {
    t <- v[[i]]$time
    o <- match(t, sort(unique(time))) # using ORDER for density plot

    d <- density(o, adjust=1) # kernel density estimation - base R function

    y <- d$y
    x <- d$x
    dd <- diff(y, lag=2)
    dt <- diff(x, lag=2)
    firstdv <- dd/dt

    n <- length(firstdv)-1
    sign <- unlist(lapply(1:n, function(i) prod(firstdv[i:i+1])))
    crit <- lapply(1:(n-1), function(i) { ifelse(prod(sign[i], sign[i+1]) < 0, i, 0)})
    cut <- which(unlist(crit)!=0)

    # Filtering to prevent noisy peaks
    if((length(cut) / 2) > length(unique(time))) {
      cut <- cut[d$y[cut] > 1/(10**10)]
    }

    ## CLUSTER DIVISION STEP (after peak detection)
    if(length(cut)==1) { # keeping cluster SAME

      new_medoid <- as.numeric(names(v)[i])
      t_groups <- v[[i]]$id
      prev_peaks[[i]] <- t_groups
      names(prev_peaks)[i] <- new_medoid

    } else {

      x_cut <- cut[c(FALSE,TRUE)] # find local min
      ncycles <- length(x_cut) + 1

      t_cut <- d$x[x_cut]
      peak_loc <- d$x[cut[c(TRUE, FALSE)]] # used to calculate period

      sep <- round(diff(peak_loc))
      sep <- which(sep >= period) # 1 check for cyclic/double-cyclic data

      # store all detected peaks in case over-clustered
      if(length(peak_loc) > 1) {
        peak_comb <- lapply(1:(length(peak_loc)-1),
                            function(i) peak_loc[c(i,i+1)])
        peak_storage[[i]] <- do.call(rbind, peak_comb) # moved from isempty(sep) to here
      }

      if(isempty(sep)) { # keeping cluster SAME
        peak_medoid <- as.numeric(names(v)[i])
        prev_peaks[[i]] <- v[[i]]$id
        names(prev_peaks)[i] <- peak_medoid
        # peak_storage[[i]] <- do.call(rbind, peak_comb) # peak storage

      } else { # divide into NEW clusters!

        breaks <- c(min(d$x), t_cut, max(d$x))
        t_groups <- lapply(1:ncycles, function(j) v[[i]]$id[which(o >= breaks[j] & o < breaks[j+1])])

        # should not contain empty groups (from noise)
        t_groups <- t_groups[!unlist(lapply(t_groups, isempty))]

        # prevent groups containing JUST 1 cell (at least 2)
        check_length <- lapply(1:length(t_groups), function(n) {
          l <- if(length(t_groups[[n]])==1 & !isempty(t_groups[[n]])) {n}
          return(l)
        })
        to_merge <- unlist(check_length)

        # determine which non-empty groups and ADD TO CLOSEST GROUP
        merge_candidates <- which(unlist(lapply(check_length, isempty)))
        merging_group <- lapply(to_merge, function(x) order(abs(merge_candidates - x))[1])

        merging_group <- unlist(merging_group)
        t_groups_0 <- t_groups
        for(clus in to_merge) {
          j <- which(to_merge %in% clus)

          t_groups[[merging_group[j]]] <- c(t_groups[[merging_group[j]]], t_groups_0[[clus]])
          t_groups[[clus]] <- NULL
        }
        t_groups <- t_groups[lapply(t_groups, length)!=1]

        peak_medoids <- lapply(t_groups, function(x) {

          points <- coordinates[x,1:2]
          a <- colMeans(points) # average point (closest to av --> new medoid)
          new_medoid <- which.min(apply(points, 1, function(p) norm(c(p-a), type='2')))
          new_medoid <- x[new_medoid]

          return(new_medoid)
        })
        names(t_groups) <- unlist(peak_medoids)
        new_peaks[[i]] <- t_groups
        cuts[[i]] <- round(abs(peak_loc))
      }
    }
  }
  all_peaks <- c(prev_peaks, unlist(new_peaks, recursive=FALSE))
  all_peaks <- all_peaks[lengths(all_peaks) != 0] # remove null entries

  cuts <- cuts[lengths(cuts) != 0]
  num_peaks <- lengths(cuts) # should have peaks corresponding to # new clusters

  # time-points that are separated by density plot
  cycle_range <- lapply(cuts, function(x) {
    if(length(x)!=2) {
      index <- cbind(1:(length(x)-1), 2:length(x))
      x <- t(apply(index, 1, function(n) x[n]))
    } else {
      x <- rbind(x, 0)
    }
    return(x)
  })

  # if a=1, then trying to separate one cycle
  ncycles <- unlist(lapply(cycle_range, function(x)
    if(sum(x[nrow(x),]) == 0) {
      a <- 1 }
    else {
      a <- nrow(x)
    }))

  # find average time-point in each cycle, use to determine period
  sep_cycles <- lapply(sort(unique(ncycles)), function(i) {
    cr <- cycle_range[ncycles==i]
    sc <- lapply(cr, function(x) x[i,])
    sc <- do.call(rbind, sc)
    sc <- matrix(sc[rowMeans(sc)!=0], ncol=2)
    return(sc) })
  sep_cycles <- sep_cycles[lengths(sep_cycles)!=0]

  if(!isempty(sep_cycles)) { # if cycles are detected post-period

    # determine dt between average time-points in respective cycles
    cycle_dt <- unlist(lapply(1:length(sep_cycles), function(i)
      mean(abs(apply(sep_cycles[[i]], 1, diff)))))

    cycle_dt <- min(ceiling(cycle_dt)) / 2 # ceiling to be conservative

  } else if(!isempty(unlist(peak_storage))) { # no cycles > period, but separation needed

    consec_dt <- apply(do.call(rbind, peak_storage), 1, diff)
    consec_dt <- round(mean(consec_dt)) + 1 # +1 to be lenient
    cycle_dt <- consec_dt

  } else { # no separation needed
    period <- NULL
    cycle_dt <- NULL
  }

  clusterLabels <- 1:nrow(coordinates)
  for (i in 1:length(all_peaks)) {
    clusterLabels[all_peaks[[i]]] <- as.numeric(names(all_peaks)[i])
  }
  return(list(clusterLabels, cycle_dt, period))
}

