#' The max-T correction
#'
#' @description Compute the max-T correction given a matrix a permuted statistics.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#'
#' @export
#' @family multcomp
compute_maxT <- function(distribution, alternative) {
  alternative <- match.arg(alternative, c("two.sided","greater","less"))
  switch (alternative,
    "greater" = {
      maxt = apply(distribution,1,max)
      p_corrected <- sapply(distribution[1,],
                            function(stati)compute_pvalue(distribution = maxt,
                                                              stat = stati,alternative = alternative))},
    "less" = {
      maxt = apply(distribution,1,min)
      p_corrected <- sapply(distribution[1,],
                            function(stati)compute_pvalue(distribution = maxt,
                                                          stat = stati,alternative = alternative))
    },
    "two.sided" = {
      maxt = apply(abs(distribution),1,max)
      p_corrected <- sapply(abs(distribution[1,]),
                            function(stati)compute_pvalue(distribution = maxt,
                                                          stat = stati,alternative = alternative))}
  )

  out <- list(main = cbind(statistic = distribution[1, ], pvalue = p_corrected))
  return(out)
}