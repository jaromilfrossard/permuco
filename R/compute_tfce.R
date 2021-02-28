#' Threshold-Free Cluster-Enhancement correction
#'
#' @description Compute the TFCE correction given a matrix a permuted statistical signals.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @param E A scalar that represent the extend parameter of the TFCE transformation. Default is \code{E = 0.5}.
#' @param H A scalar that represent the height parameter of the TFCE transformation. Default is \code{H = 1}.
#' @param ndh The number of terms in the approximation of the integral.
#' @export
#' @family multcomp
compute_tfce <- function(distribution, alternative = "greater", E = 0.5, H = 1, ndh){
  alternative <- match.arg(alternative, c("two.sided","greater","less"))


  range_d <- range(abs(distribution))
  dhi <- seq(from = range_d[1], to =range_d[2], length.out = ndh)
  dh <- dhi[2]-dhi[1]

  # distribution_tfce <- apply(abs(distribution),1, function(rowi){
  #   max(compute_tfce_statistic(signal = rowi, E = E, H = H,dh = dh, dhi = dhi))})

  if(alternative=="two.sided"){
    distribution_tfce <- tfce_distribution(abs(distribution), E, H, dh, dhi)
  }else if(alternative=="greater"){
    distribution_tfce <- tfce_distribution(distribution, E, H, dh, dhi)
  }else if(alternative=="greater"){
  distribution_tfce <- tfce_distribution(-distribution, E, H, dh, dhi)
  }
mx_tfce <-apply(distribution_tfce,1,max)



  #statistic <- compute_tfce_statistic(signal = abs(distribution[1,]), E = E, H = H,dh = dh, dhi = dhi)
  statistic <-   distribution_tfce[1,]


  pvalue <- sapply(statistic,function(s)compute_pvalue(distribution = mx_tfce,stat=s))

  out <- list(main = cbind(statistic = statistic*sign(distribution[1,]), pvalue = pvalue),
              distribution = mx_tfce, E = E, H = H,
              dhi = dhi)
  return(out)



}


