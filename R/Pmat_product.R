#' Multiplies a vector with a Pmat object
#'
#' @description Multiplies a vector or matrix with a Pmat object
#' @param x A a vector or matrix
#' @param P A Pmat object
#' @param type A character string indicating the type of transformation. This argument need to be specified if P is not of class Pmat.
#' @return A matrix n x np containing the permutated/signflipped vectors.
#' @export
Pmat_product <- function(x, P, type = NULL){UseMethod("Pmat_product")}

#' @export
Pmat_product.numeric <- function(x, P, type = NULL){
  if(is.null(type)){type = attr(P,"type")}
  ### check length
  if(length(x)!=nrow(P)){warning("x and P should have matching row length.")}
  ### product
  switch(type,
         "permutation" ={x <- matrix(x[as.matrix(P)],ncol = np(P))},
         "signflip" = {x <- matrix(x*as.matrix(P),ncol = np(P))}
         )
  return(x)
}

#' @export
Pmat_product.matrix <- function(x, P, type = NULL){
  if(is.null(type)){type = attr(P,"type")}
  ## check x is matrix
  if(ncol(x)==1){return(Pmat_product(x = as.numeric(x), P = as.matrix(P), type = type))}
  ### check P is vector
  if(is.matrix(P)){
    if(ncol(P)!=1){warning("Only 1 transformation is allowed when x is a matrix.")}
    else{
      P <- as.numeric(P)}}
  ### check length
  if(nrow(x)!=length(P)){warning("x and P should have matching row length.")}
  ## product
  switch(type,
         "permutation" ={x <- matrix(x[as.numeric(P),],ncol = ncol(x))},
         "signflip" = {x <- matrix(x*as.numeric(P),ncol = ncol(x))})

  return(x)
}
