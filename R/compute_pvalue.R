compute_pvalue <- function(distribution, stat = distribution[1], alternative="two.sided", na.rm = T, digits = 10){
  if(na.rm){distribution = as.numeric(na.omit(distribution))}
  distribution=round(distribution,digits = digits)
  stat=round(stat, digits = digits)
  stat=matrix(stat,nrow = 1)
  switch(alternative,
         "two.sided" = {apply(stat,2, function(val)mean(abs(distribution) >= abs(val) , na.rm = na.rm))},
         "less" = {apply(stat,2, function(val)mean(distribution <= val , na.rm = na.rm))},
         "greater" = {apply(stat,2, function(val)mean(distribution >= val, na.rm = na.rm))})
}