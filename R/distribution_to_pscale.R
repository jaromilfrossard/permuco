distribution_to_pscale <- function(distribution, test, alternative){
  if(test == "fisher"){
    out <- apply(distribution,2,function(col){compute_all_pvalue(col,alternative = "less")})
  }else if (test == "t"){
    switch(alternative,
           "greater" = {out <- apply(distribution,2,function(col){compute_all_pvalue(col,alternative = "less")})},
           "less" = {out <- apply(distribution,2,function(col){compute_all_pvalue(col,alternative = "greater")})},
           "two.sided" = {out <- apply(abs(distribution),2,function(col){compute_all_pvalue(col,alternative = "less")})})}
  return(out)
}