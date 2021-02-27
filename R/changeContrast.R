#'@importFrom stats contrasts<-
changeContrast <- function(data,contr){
  isfact <- sapply(data, function(c){is.factor(c)||is.character(c)})
  for(i in c( 1:ncol(data))[isfact]){
    data[ ,i] <- as.factor(data[ ,i])
    data[ ,i] <- droplevels(data[ ,i])
    contrasts(data[ ,i]) <- contr
  }
  return(data)
}