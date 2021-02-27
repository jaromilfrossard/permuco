#' Method to convert into \code{Pmat} object.
#'
#'@description Convert a matrix into a \code{Pmat} object.
#'
#'@param x a matrix.
#'
#'@export
#'@family pmat
as.Pmat <- function(x) {UseMethod("as.Pmat")}

#' @export
#' @family pmat
as.Pmat.matrix <- function(x){
  np = NCOL(x)
  n=NROW(x)
  v=1:n
  #check thirst column
  if(sum(x[, 1] == v) != n){stop("cannot be coherce into a Pmat object : the first row should be a 1:n vector")}
  #check the rest
  if(sum(apply(x[, -1],2 ,function(p){sum(sort(p) == v) == n})) != np-1){
    stop("cannot be coherce into a Pmat object : the matrix should be compose of permutation of the 1:n vector")
  }
  attr(x, which = "type") = "default"
  attr(x, which = "np") = np
  class(x) <- "Pmat"
  return(x)
}

#' @export
#' @family pmat
as.matrix.Pmat <- function(x, ...){
  return(matrix(x, ncol = NCOL(x)))
}
