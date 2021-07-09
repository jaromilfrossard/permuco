############################################################################################################
#' Clustermass test correction
#'
#' @description Compute the clustermass test correction given a matrix a permuted statistical signals.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param threshold A scalar that respresents the threshold to create the clusters.
#' @param aggr_FUN A function to compute the clustermasses. See details for examples.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#'
#' @details The \code{aggr_FUN} argument may take predefined function as the sum: \code{aggr_FUN = sum} and also user-defined function as the sum of squares: \code{aggr_FUN = function(x){sum(x^2)}}
#' @export
#' @family multcomp
compute_clustermass <-  function(distribution, threshold, aggr_FUN, alternative = "greater"){
  switch(alternative,
         "greater" = {
           #threshold <- abs(threshold)
           #selected <- distribution > threshold
           extreme <- function(x)max(x,na.rm = T)
         },
         "less" = {
           #threshold <- -abs(threshold)
           #selected <- distribution < threshold
           extreme <- function(x)max(x,na.rm = T)},
         "two.sided" = {
           distribution <- abs(distribution)
           # threshold <- abs(threshold)
           #selected <- distribution > threshold
           extreme <- function(x)max(x,na.rm = T)
         })

  ##create connected component labeling
  # cl <- (selected-cbind(0,selected[,-NCOL(selected),drop=F]))==1
  # cl <- t(apply(cl,1,cumsum))
  # cl[!selected] <- 0
  cl <- get_cluster(distribution = distribution, threshold = threshold, alternative = alternative, side= "all")

  mass_distribution <- sapply(1:NROW(cl),function(rowi){
    extreme(sapply(1:max(1,max(cl[rowi,])),function(i){
      aggr_FUN(distribution[rowi,cl[rowi,] == i])
    }))
  })
  mass_statistic <- sapply(1:max(cl[1,]),function(i){
    aggr_FUN(distribution[1,cl[1,] == i])
  })
  statistic <- c(NA,mass_statistic)[cl[1,]+1]


  pvalue <- sapply(mass_statistic,function(mi)
    compute_pvalue(stat = mi,distribution = mass_distribution, alternative = "greater"))
  pvalue <- c(NA,pvalue)[cl[1,]+1]
  out <- list(main = cbind(statistic = statistic, pvalue = pvalue, cluster_id = cl[1,]),
              distribution = mass_distribution,
              threshold = threshold)
  return(out)
}