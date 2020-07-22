#' @importFrom stats na.omit
compute_all_pvalue <- function(distribution, alternative="two.sided", na.rm = T, digits = 10){
  if(na.rm){distribution = na.omit(distribution)}
  distribution=round(distribution,digits = digits)
  np = length(distribution)
  switch(alternative,
         "two.sided" = {distribution = - abs(distribution)},
         "less" = {},
         "greater" = {distribution = - distribution})
  ceiling(rank(distribution))/np
}