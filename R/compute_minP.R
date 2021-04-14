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

  minimal_pval <- table(minp)[1]/nrow(distribution)
  if(minimal_pval>0.02)warnings(paste0("The minimal adjusted p-value is ",minimal_pval,". Try to increase the number of permutations or you compare to much uncorrelated variables."))

  pvalue <- distribution_rank[1,]
  p_corrected <- sapply(pvalue,function(pi)compute_pvalue(distribution = minp,
                                                          stat = pi,alternative = "less"))
  out <- list(main = cbind(statistic = distribution[1, ], pvalue = p_corrected))
  return(out)
}