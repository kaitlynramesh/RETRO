# Pseudotime fitting function
pseudotime_fit <- function(retro_pt_obj, retro_obj, h=0.5) {
  # Obtain coordinates and lineage information 
  coordinates <- retro_pt_obj$Coordinates
  time <- retro_pt_obj$Time
  cells_to_lin <- retro_pt_obj$Cells_to_Lin
  lin_membership <- retro_pt_obj$Lin_Membership
  nl <- retro_pt_obj$Num_Lineages
  
  # Extract MST specifics
  cell_MST <- retro_obj$RETRO_MST
  id <- cell_MST$ID
  lineages <- cell_MST$Lineages
  clusterLabels <- cell_MST$ClusterLabels
  reclus <- cell_MST$Reclus

  # Initialize lists for storage
  pseudotime_list <- vector(mode='list', length=nl)
  lambda_list <- vector(mode='list', length=nl)
  arclen_list <- vector(mode='list', length=nl)
  graph_list <- vector(mode='list', length=nl)
  
  # Obtain curve and arc length information
  ext_bcurves <- sapply(retro_pt_obj$Curve, "[", 1)
  arclengths <- sapply(retro_pt_obj$Curve, "[", 4)
  
  # Initialize lists for storage
  pseudotime_list <- vector(mode='list', length=nl)
  arclen_list <- vector(mode='list', length=nl)
  graph_list <- vector(mode='list', length=nl)
  
  # Obtain curve and arc length information
  ext_bcurves <- sapply(retro_pt_obj$Curve, "[", 1)
  arclengths <- sapply(retro_pt_obj$Curve, "[", 4)
  
  # Project cells according to the corresponding MST segment 
  projection_res = projection_by_segment(nl=nl, 
                                         coordinates=coordinates,
                                         lineages=lineages,
                                         lin_membership=lin_membership,
                                         clusterLabels=clusterLabels,
                                         id=id, arclengths=arclengths,
                                         ext_bcurves=ext_bcurves)
  
  lambda_per_lineage = sapply(projection_res, "[", "final_lambda")
  time_per_lineage = sapply(projection_res, "[", "real_time")
  which_cells = sapply(projection_res, "[", "which_cells")
  
  # Pseudotime-fitting per lineage (for cyclic/multifurcating lineages)
  for(i in 1:nl) {
    real_time = time_per_lineage[[i]]
    lambda = lambda_per_lineage[[i]]
    
    curve = cbind(Arc_Length=lambda, Time=real_time)
    curve_fit = lambda_curve_fit(lambda=lambda, 
                                 real_time=real_time)
    
    wc = which_cells[[i]][order(real_time)]
    pseudotime_list[[i]] = populate_matrix(nrow(coordinates),
                                           wc, 
                                           curve_fit[,2])
    
    lambda_list[[i]] = populate_matrix(nrow(coordinates),
                                       wc, 
                                       curve_fit[,1])
    
    g <- ggplot(as.data.frame(curve), aes(x=Arc_Length, y=Time)) + 
      geom_point() +
      geom_point(data=curve_fit, aes(x=x,y=y), col='lightgreen') +
      theme_bw() + 
      theme(axis.text=element_text(size=10))
    graph_list[[i]] <- g
    
  }

  # Average shared pseudotime values
  if(length(lineages)>1) {
    pseudotime_mat <- do.call(cbind, pseudotime_list)
    pseudotime <- rowMeans(pseudotime_mat, na.rm=TRUE)
  } else {
    pseudotime_mat = as.numeric(pseudotime_list[[1]])
    pseudotime <- pseudotime_list[[1]]
  }
  
  retro_pt_obj$Pseudotime <- pseudotime
  retro_pt_obj$Pseudotime_Mat = pseudotime_mat
  retro_pt_obj$Graph_List <- graph_list
  retro_pt_obj$Fitting <- lambda_list
  
  return(retro_pt_obj)
}


lambda_curve_fit <- function(lambda, real_time) {
  
  med_fitting = as.data.frame(cbind(lambda, real_time))
  colnames(med_fitting) = c("x", "y")
  lambda_by_time = split(lambda, real_time)
  t = as.numeric(names(lambda_by_time))
  
  med_lambda = lapply(lambda_by_time, function(x) as.numeric(quantile(x)[2:4]))
  med_lambda = cbind(unlist(med_lambda), sort(rep(t, 3)))
  
  colnames(med_lambda) = c("x", "y")
  f = monpol(y ~ x, data=as.data.frame(med_lambda), deg=3)
  predicted_values <- predict(f, newdata = med_fitting)
  
  curve = data.frame(cbind(med_fitting[,1], predicted_values))
  colnames(curve) = c("x", "y")
  
  return(curve)
}

scale_time <- function(time, pseudotime) {
  or <- as.vector(range(pseudotime)) # original range
  ar <- as.vector(range(time)) # actual range
  scaled_pt <- lapply(pseudotime, function(p) ((((p-or[1])*(ar[2]-ar[1])) / (or[2]-or[1])) + ar[1]) )
  scaled_pt <- unlist(scaled_pt)
  return(scaled_pt)
}

# Create matrix of pseudotime values (1) if cell has >1 PT value and (2) to ensure that
# correct pseudotime assigned to each cell
populate_matrix <- function(max_cells, l, pt) {
  pseudotime_mat <- as.matrix(rep(NA, length=max_cells))
  for (i in 1:length(l)) {
    pseudotime_mat[l[i],] <- pt[i]
  }
  return(pseudotime_mat)
}

# Update lineage information such that EVERY cell maps to a lineage
get_mapped_cells <- function(retro_obj) {
  
  coordinates <- retro_obj$Coordinates
  time <- retro_obj$Time
  lin_membership <- retro_obj$Lin_Membership
  cells_to_lin <- retro_obj$Cells_to_Lin
  cell_MST <- retro_obj$RETRO_MST
  clusterLabels <- cell_MST$ClusterLabels # for node to lineage consistency
  lineages <- cell_MST$Lineages
  id = cell_MST$ID
  
  # Identify the cells that don't belong to a lineage
  total_cells <- 1:nrow(coordinates)
  x <- total_cells[!(total_cells %in% unique(unlist(lin_membership)))]
  nl <- length(cells_to_lin)
  
  if (length(x)==0) {
    
    final_cells_to_lin <- cells_to_lin
    final_lin_mem <- lin_membership
    
  } else if (nl==1) {
    
    final_cells_to_lin = list(coordinates)
    final_lin_mem = list(1:nrow(coordinates))
    
  } else {
    
    unmatched_clusters = unique(clusterLabels[x])
    final_lin_mem = lin_membership # initialize list
    
    cluster_membership <- lapply(lin_membership, function(l) as.numeric(names(table(clusterLabels[l]))))
    lin_ids = unique(unlist(cluster_membership))
    
    nearest_node = lapply(unmatched_clusters, function(u) {
      
      p = coordinates[u,]
      dist_mat = as.matrix(dist(rbind(p, coordinates[lin_ids,])))
      dist_mat = dist_mat[-1,]
      dist_vec = dist_mat[,1]
      
      return(lin_ids[which.min(dist_vec)])
    })
    nearest_node = as.numeric(nearest_node)
    
    nearest_lin <- sapply(nearest_node, function(x) 
      which(sapply(cluster_membership, function(y) x %in% y)))
    
    cluster_match_df = cbind(unmatched_clusters, nearest_node)
    
    for(r in 1:nrow(cluster_match_df)) {
      u = cluster_match_df[r,1]
      l = nearest_lin[[r]]
      to_match = x[which(clusterLabels[x] %in% u)]
      final_lin_mem[l] = lapply(final_lin_mem[l], function(x) c(x, to_match))
    }
    
    final_cells_to_lin <- lapply(final_lin_mem, function(x) coordinates[x,])
  } 
  
  # Initialize pseudotime object 
  retro_pt_obj <- retro_obj
  retro_pt_obj$Lin_Membership <- final_lin_mem
  retro_pt_obj$Cells_to_Lin <- final_cells_to_lin
  retro_pt_obj$Graph_List <- NULL
  retro_pt_obj$Arc_Length <- NULL
  retro_pt_obj$Pseudotime <- NULL
  
  return(retro_pt_obj)
}


projection_by_segment <- function(nl, coordinates, lineages, lin_membership, clusterLabels, id, arclengths, ext_bcurves) {
  res_all = vector(mode="list", length=nl) # initialize list for projection data 
  for(i in 1:nl) {
    # Obtain cluster/nodes belonging to each lineage
    l <- unlist(lineages[i])
    lin_clusterLabels = clusterLabels[lin_membership[[i]]] 
    lin_id = id[l]
    
    # Group nodes together if not in lineage
    node_not_used = unique(lin_clusterLabels[!(lin_clusterLabels %in% id[l])]) # not in lineage
    node_labels = lin_clusterLabels
    
    # f cluster label is not a part of lineage, map to existing node in lineage
    if(!isempty(node_not_used)) {
      nodes_to_match = lin_id[-which(lin_id %in% node_not_used)] # nodes in final lineage
      if(isempty(nodes_to_match)) { nodes_to_match = lin_id }
      
      node_match = sapply(node_not_used, function(node) {
        coord1 = cbind(coordinates[,1:2], time)[node,] 
        coord2 = cbind(coordinates[,1:2], time)[nodes_to_match,] 
        
        dist_mat = as.matrix(dist(rbind(coord1, coord2)))
        dist_mat = dist_mat[-1,]
        closest_node = nodes_to_match[which.min(dist_mat[,1])] # closest node in lineage
        return(closest_node)
      }) 
      for(n in seq(length(node_match))) { # replace cluster label with matched node
        node_labels[node_labels==node_not_used[n]] = as.numeric(node_match[n])
      }    
    }

    # Obtain list of cells belonging to each node of MST lineage
    node_ids = lapply(split(lin_membership[[i]], node_labels), as.numeric)
    num_nodes = length(node_ids)
    
    # Edge list for order of nodes in lineage
    edge_index <- apply(cbind(l[1:(length(l)-1)], l[-1]), 2, function(x) id[x])
    
    # Dimensions used for projection
    dimensions = 1:2 # c(1:3) # c(1:2, ncol(coordinates))
    
    # Determine portion of curve corresponding to each edge
    curve_per_edge <- lapply(0:(nrow(edge_index)+1), function(j) { 
      section <- seq( 1001*j, 1001*(j+1), 1) # 1001 (=seq(0,1,.001)) points per edge
      
      curve_section <- ext_bcurves[[i]][section,dimensions] #1:2]
      return(curve_section)
    })
    
    # Node labels/coordinates in order of lineage
    nodes_in_lin <- c(edge_index[,1], tail(as.vector(edge_index), 1)) 
    node_coord <- lapply(node_ids, function(x) coordinates[x,dimensions])
    
    # Obtain projection of cells to corresponding segment of lineage 
    projection <- lapply(1:length(nodes_in_lin), function(j)  {
      node <- nodes_in_lin[j]
      ids = node_ids[[as.character(nodes_in_lin[j])]] 
      node_cluster = node_coord[[as.character(nodes_in_lin[j])]]

      if(is.null(node_cluster)) {
        ids = NULL
        scaled_lambda = NULL
      } else {
        if(j==1) { # starting node projects to extended curve
          segment <- do.call(rbind, curve_per_edge[c(j,j+1,j+2)])
        } else if (j==num_nodes) { # ending node projects to extended curve
          segment <- do.call(rbind, curve_per_edge[c(j-1,j,j+1)])
          proj <- project_to_curve(x=node_cluster, s=segment, stretch=2)
        } else { # regular nodes project to major chunks of curve
          segment <- do.call(rbind, curve_per_edge[c(j-1,j,j+1)])
        }

        node_cluster = rbind(segment[1,], node_cluster) # edit to accommodate clusters w/ 1 cell
        proj <- project_to_curve(x=node_cluster, s=segment, stretch=2)
        proj_lambda = proj$lambda[-1]
        
        # scale lambda to corresponding ARC LENGTH (prev to current)
        or <- as.vector(range(proj_lambda))
        ar <- c(sum((head(arclengths[[i]], j-1))), sum((head(arclengths[[i]], j))))
        if(diff(or)!=0) {
          ar <- c(sum((head(arclengths[[i]], j-1))), sum((head(arclengths[[i]], j))))
          scaled_lambda <- lapply(proj_lambda, function(l) ((((l-or[1])*(ar[2]-ar[1])) / (or[2]-or[1])) + ar[1]) )
          scaled_lambda <- unlist(scaled_lambda)
        } else {
          scaled_lambda <- rep(mean(ar), length=length(proj_lambda))
        }
      }
      return(list(ids, scaled_lambda))
    })
    
    which_cells <- unlist(sapply(projection, "[", 1)) 
    final_lambda <- unlist(sapply(projection, "[", 2))
    real_time <- time[which_cells]
    
    x = list("final_lambda"=as.numeric(final_lambda),
             "real_time"=real_time,
             "which_cells"=which_cells)
    res_all[[i]] = x
  }
  return(res_all)
}
