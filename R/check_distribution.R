check_distribution <- function(distribution, digits = 10, n_unique = 30){
  unique_d = apply(distribution,2,function(d){
    length(unique(round(d,digits = digits)))
  })
  if(sum(unique_d<n_unique)>0){
    warning(paste("the distribution of ",paste(colnames(distribution)[unique_d<n_unique],sep =", ",collapse = ", "), " may be discrete.", sep=" ",collapse = " "))
  }
}