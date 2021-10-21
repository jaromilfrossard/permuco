###########################

#' Get the ID of the cluster on a distribution matrix
#'
#' @description Get if of clusters given a threshold
#'
#' @param distribution a matrix of the distribution.
#' @param threshold a numeric indicating the threshold.
#' @param alternative A character string indicating the alternative hypothesis. Choose between \code{"two.sided"}, \code{"greater"}, \code{"less"}. Default is \code{"two.sided"}.
#' @return a matrix with integer indicating the id of the clusters
#' @useDynLib permuco
#' @importFrom Rcpp sourceCpp
#' @keywords internal
get_cluster <- function(distribution, threshold, alternative = "two.sided", side = "all"){
  UseMethod("get_cluster")
}

get_cluster.matrix <- function(distribution, threshold, alternative, side){
  alternative <- match.arg(alternative, c("two.sided","greater","less"))
  side <- match.arg(side, c("all","starting","ending"))


  if(alternative == "two.sided"){
    distribution <- abs(distribution)
    threshold <- abs(threshold)
}
  if(alternative == "greater"){
    distribution <- distribution
    threshold <- abs(threshold)
    }
  if(alternative == "less"){
    distribution <- -distribution
    threshold <- abs(threshold)
  }

  cl <- get_cluster_matrix(distribution, abs(threshold), side = side)

}