compute_benjamini_hochberg <- function(statistic = NULL, pvalue){
  n <- length(pvalue)
  i <- n:1L
  o <- order(pvalue, decreasing = TRUE)
  ro <- order(o)
  pvalue <- pmin(1, cummin(n/i * pvalue[o]))[ro]
  out <- list(main = cbind(statistic = statistic,pvalue = pvalue))
  return(out)
}