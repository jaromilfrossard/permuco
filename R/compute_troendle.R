#' The Troendle's correction
#'
#' @description Compute the Troendle's correction given a matrix a permuted statistics.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#'
#' @export
#' @family multcomp
compute_troendle = function (distribution, alternative) {
  distribution_rank <- apply(distribution,2,function(col){compute_all_pvalue(col,alternative = alternative)})
  length_unique_minp <- length(unique(apply(distribution_rank,1,min)))

  if(length_unique_minp<200)warnings(paste0("The adjusted distribution contains only ",length_unique_minp," unique values. Try to increase the number of permutation."))

  pvalue <- distribution_rank[1,]

  rank_uncorr <- rank(pvalue)


  # order test from low to high
  test_order <- sort(unique(rank_uncorr))
  order_test_list <- lapply(1:length(test_order),function(testi){
    ri <- test_order[testi]
    list(infos = cbind(order=testi,col=which(rank_uncorr==ri)),
         mins = apply(distribution_rank[,which(rank_uncorr==ri),drop=F],1,min))
  })

  # take cumulative minimum of the ordered distributions of p's
  cummins <- do.call("cbind",lapply(order_test_list,function(oi)(oi$mins)))
  cummins <- lapply(1:nrow(cummins),function(rowi){rev(cummin(rev(cummins[rowi,,drop=F])))})
  cummins <- do.call("rbind",cummins)

  # Compute p-value

  p_corrected <- apply(cummins,2,
                       function(coli)compute_pvalue(coli,alternative = "less"))

  # Compute inverse order

  infos <- do.call("rbind",lapply(order_test_list,function(oi)(oi$infos)))

  # Inverse order of the increasing p-values.
  p_corrected <- (cummax(p_corrected) [infos[,1]])[order(infos[,2])]

  out <- list(main = cbind(statistic = distribution[1, ], pvalue = p_corrected))
  return(out)
}