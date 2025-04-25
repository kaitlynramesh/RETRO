#' @importFrom methods setClass new
#'
#' @title RETRO: S4 class for pseudotime analysis.
#' @description This S4 class contains slots for containing data and parameters necessary
#' for analysis with RETRO functions, and also stores the final inferred pseudotime.
#' @aliases RETRO
#' @export
#'
#' @slot coordinates Coordinates determined by \code{weight_coord()} that include a contribution from sampling time.
#' @slot coordinates_uw Coordinates determined by \code{weight_coord()} that have an **optional** contribution from sampling time.
#' It is recommended that these coordinates are un-weighted so that the initial clustering covers the entire low-dimensional projection.
#' @slot time Numeric vector of sampling time points corresponding to each cell.
#' @slot k_range Numeric vector for the range of K clusters that will be used to partition the coordinates
#' @slot start (optional) Character that is either "Average" or "Mode" to determine the starting node (if no starting cells are specified).
#' If "Average," the starting node is assigned to the cluster with the lowest average sampling time.
#' If "Mode," the starting node is assigned the cluster with the greatest representation of the lowest sampling time-point.
#' Default is "Mode".
#' @slot terminal_cells (optional) List of cells belonging to the terminal ends of the trajectory.
#' Each list element can contain the cells corresponding to one terminal cell type.
#' @slot starting_cells (optional) Vector of cells belonging to the starting end of the trajectory.
#' RETRO currently allows for only one starting node.
#' @slot threshold (optional) Cutoff to classify whether clusters are terminal ends of the trajectory. Default is 0.1.
#' Note that this cutoff is scaled by number of clusters (K) / maximum number of clusters specified (Max_K).
#' @slot period (optional) Specifies the minimum difference in time where one can expect the cells
#' to return to an earlier state gene expression (make a cycle).
#'
#' @slot all_k List containing cluster labels and centroids for each iteration of clustering on the range of K
#' @slot all_scores List containing MST scores and cell membership for all K clusters in the specified range
#' @slot num_lineages !!! include description here
#' @slot centroids !!!! include description here
#' @slot lin_membership !!!!!!!description goes here
#' @slot cells_to_lin !!!!!!description goes here
#' @slot RETRO_MST !!!!!!!description goes here
#' @slot RETRO_Curve !!!!!!description goes here
#' @slot pseudotime_graph !!! include description here
#' @slot pseudotime !!!! include description here
#'
#' @return An object of class \code{RETRO}
#'
#'
set_RETRO_class <- setClass("RETRO", slots=c(coordinates="matrix",
                                             coordinates_uw="matrix",
                                             time="numeric",
                                             k_range = "numeric",
                                             all_k = "list",
                                             all_scores = "list",
                                             period="ANY",
                                             threshold = "numeric",
                                             starting_cells = "ANY",
                                             terminal_cells= "ANY",
                                             start = "character",
                                             num_lineages = "numeric",
                                             centroids = "list",
                                             lin_membership = "list",
                                             cells_to_lin = "list",
                                             RETRO_MST = "list",
                                             RETRO_Curve = "list",
                                             pseudotime_graph = "list",
                                             pseudotime = "numeric"),
                            prototype=list(threshold=0.1,
                                           start = "Mode",
                                           starting_cells=NULL,
                                           terminal_cells=NULL,
                                           period=NULL))


