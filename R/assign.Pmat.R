
#' @export
"[<-.Pmat" <- function(x,i,j,value){
  x <- as.matrix(x)
  x[i,j] <- value
  return(x)
}


#' @export
`[.Pmat` <- function(x,i,j,drop = FALSE){
  if(drop){warning("drop is set to TRUE")}
  mc <- match.call()
  attr <- attributes(x)
  x <- as.matrix(x)
  x <- x[i,j,drop = FALSE]
  if(is.null(mc$i)){
    attr(x,"type") <- attr$type
    attr(x,"counting") <- attr$counting
    attr(x,"np") <- ncol(x)
    attr(x,"n") <- attr$n
    attr(x,"class") <- attr$class}
  return(x)
}