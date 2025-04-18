library(bezier)
library(MonoPoly)
library(ggplot2)
library(princurve)
library(RColorBrewer)
library(scales)
library(viridis)

### POST CLUSTER ###
create_meta_obj <- function(obj, coordinates, terminal_cells=NULL, starting_cells=NULL, 
                            start=NULL, period=NULL, k,
                            clustering, scoring) {
  time <- obj$Time
  meta_obj <- list('Coordinates'=coordinates,
                   'Time'=time, 'Period'=period,
                   'Terminal_Cells'=terminal_cells, 'Starting_Cells'=starting_cells, 'Start'=start,
                   'K'=k, 'All_K'=clustering, 'All_Scores'=scoring,
                   'Num_Lineages'=NULL, 
                   'Lin_Membership'=NULL, 'Cells_to_Lin'=NULL,
                   'RETRO_MST'=NULL, 'Curve'=NULL)
  return(meta_obj)
}


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



# Determines # of Lineages (Heatmap) and Cell Membership to Each Lineage
get_num_lineages <- function(retro_obj, percent = 0.05, cutoff=0.8, threshold) {
  
  all_k <- retro_obj$All_K
  all_scores <- retro_obj$All_Scores
  
  coordinates <- retro_obj$Coordinates
  time <- retro_obj$Time
  start <- retro_obj$Start
  period <- retro_obj$Period
  starting_cells <- retro_obj$Starting_Cells
  terminal_cells <- retro_obj$Terminal_Cells
  max_k <- max(retro_obj$K)
  num_scores <- length(all_scores[[1]])
  
  top_data <- get_top_scores(all_scores, num_scores, percent)
  membership <- top_data[[1]] # membership matrices
  top_scores = top_data[[2]]
  circular <-  isempty(unlist(lapply(membership, nrow))) # must all be ==1 
  
  nl_detected <- unlist(lapply(membership, function(i) {  
    x <- nrow(i)
    nl <- ifelse(is.null(x), 1, x)
    return(nl)
  })) # obtain lineages detected from top MSTs

  # (1) CLUSTER MEMBERSHIP MATRIX
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
  
  # (2) CLUSTER CELLS BASED ON MEMBERSHIP MATRIX
  factor = pam(m, nl_knn, diss=TRUE) # grouping for lineage membership matrices
  lin_membership <- lapply(1:nl_knn, function(i) { 
    x <- t(all_membership[,factor$clustering==i])
    return(which(colMeans(x) >= cutoff))})
  
  # (3) IDENTIFY CELLS THAT MAP TO LINEAGE
  binary_mat = matrix(0, nrow=nrow(coordinates), ncol=length(lin_membership))
  for (j in seq_along(lin_membership)) {
    binary_mat[lin_membership[[j]], j] <- 1
  }
  
  # (4) IDENTIFY CLOSEST MST TRAJECTORY
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
  
  index = which.min(closest_mst) # based on membership matrix
  
  f <- top_scores[index]
  clusn.i <- ifelse (f %% num_scores == 0, f / num_scores, floor(f / num_scores) + 1) # which k
  iteration.i <- ifelse (f %% num_scores == 0, num_scores, f %% num_scores) # which iterations
  
  cell_MST <- create_dMST(coordinates, all_k[[clusn.i]][[iteration.i]], time, terminal_cells, 
                          starting_cells, threshold, max_k=max_k, start=start, period=period)
  id <- cell_MST$ID 
  nl <- length(cell_MST$Lineages) # final # of lineages
  
  centers <- lapply(1:nl, function(i) coordinates[id[unlist(cell_MST$Lineages[[i]])],])
  
  retro_obj$Num_Lineages <- nl
  retro_obj$RETRO_MST <- cell_MST
  retro_obj$Centroids <- centers
  
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
  
  retro_obj$Cells_to_Lin <- cells_to_lin
  retro_obj$Lin_Membership <- lin_membership
  retro_obj$Top_Data <- top_data
  
  return(retro_obj)
}




### BEZIER CURVE CONTROL POINTS ###############################################

# Extension for Bézier curve to cover all points
extend_centers <- function(centroids, p) {
  # starting center extension
  v_start <- centroids[1,] - centroids[2,]
  root_center <- centroids[1,] + p * v_start
  
  # ending center extension
  n <- nrow(centroids)
  v_terminal <- centroids[n,] - centroids[n-1,]
  terminal_center <- centroids[n,] + p * v_terminal
  
  extreme_centers <- rbind(root_center, terminal_center)
  
  return(extreme_centers)
}

# Contains concatenated Bézier curve functions, manages extension of Bézier curve 
get_bezier_curve <- function(retro_obj, extension=2) {
  
  nl <- retro_obj$Num_Lineages
  
  bcurve_data <- lapply(1:nl, function(i) {
    centroids <- retro_obj$Centroids[[i]]
    
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
  
  retro_obj$Curve <- bcurve_data
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


#### VISUALIZATION ############################################################
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

### K-Medoids Scoring Box Plots
boxplot_scoring <- function(all_scores, num_scores, num_k, k) {
  
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
    # geom_text(aes(label=round(av_time, 3)), size=2.5, hjust=1.5, check_overlap=TRUE, color='black') +
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


# Benchmark w/ Psupertime Density Graph 
psupertime_density <- function(x, time, pseudotime, bw=0.5) {
  
  # x: coordinates of low-dim projection
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
  
  # Step 1: Create a table of custom bandwidths
  bw_list <- proj_dt %>%
    group_by(label_input) %>%
    summarize(bw = bw, .groups = "drop")  
  
  # Step 2: Nest data and compute density per group using your bw_list
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
  print(bw_list_new)
  
  
  ggplot(dens_df, aes(x = x, y = y, color = label_input, fill = label_input)) +
    geom_line() +
    geom_area(alpha = 0.4, position = "identity") +
    scale_colour_manual(values = col_vals) +
    scale_fill_manual(values = col_vals) +
    # scale_y_continuous(trans = 'log1p') +
    labs(
      x = "Pseudotime",
      y = "Density",
      color = "Time",
      fill = "Time"
    )

}
