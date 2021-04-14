#' The min-P correction
#'
#' @description Compute the min-P correction given a matrix a permuted statistics.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#'
#' @export
#' @family multcomp
compute_minP <- function(distribution, alternative) {
  distribution_rank <- apply(distribution,2,function(col){compute_all_pvalue(col,alternative = alternative)})
  minp <- apply(distribution_rank,1,min)

  length_unique_minp <- length(unique(minp))
  if(length_unique_minp<200)warnings(paste0("The adjusted distribution contains only ",length_unique_minp," unique values. Try to increase the number of permutation."))

  pvalue <- distribution_rank[1,]
  p_corrected <- sapply(pvalue,function(pi)compute_pvalue(distribution = minp,
                                                          stat = pi,alternative = "less"))
  out <- list(main = cbind(statistic = distribution[1, ], pvalue = p_corrected))
  return(out)
}