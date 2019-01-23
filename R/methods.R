#PRINT methods=============================

#' @export
print.lmperm <- function(x,...){
  print(x$table,...)
}

#' @export
print.lmpermutation_table<-function(x, digits = 4, na.print = "", ...){
  cat(attr(x,"heading"), sep = "\n\n")
  cat(attr(x,"type"), sep = "\n\n")
  if(is.data.frame(x)){
    print(as.matrix(x), digits = digits, na.print = na.print, ...)
  }else{
    print.listof(x, digits = digits, ...)
  }
}

#PLOT methods=============================

#'Plot method for class \code{"lmperm"}.
#'
#'@description Show the density of statistics and the test statistic.
#'
#'@param x A \code{"lmperm"} object.
#'@param FUN A function to compute the density. Default is \link{density}.
#'@param ... futher arguments pass to plot.
#'@details Other argument can be pass to the function : \cr \cr
#'\code{effect} : a vector of character string indicating the name of the effect to plot. \cr
#'@importFrom graphics par plot abline
#'@importFrom stats density
#'@export
plot.lmperm <- function(x, FUN = density, ...){
  par0 = par()
  #data
  distr = x$distribution

  dotargs = list(...)
  if(is.null(dotargs$effect)){
    effect <- colnames(distr)
  }else{effect <- colnames(distr)[which(colnames(distr)%in%dotargs$effect)]}
  distr = distr[,which(colnames(distr)%in%effect),drop=F]
  dotargs = dotargs[names(dotargs)!="effect"]

  #subplot
  p <- NCOL(distr)
  div <- seq_len(abs(p))
  factors <- div[p %% div == 0L]
  mfrow1 <- factors[ceiling(length(factors)/2)]
  mfrow <- c(mfrow1,p/mfrow1)
  par(mfrow = mfrow)

  #plot
  for(i in 1:NCOL(distr)){
    plot(FUN(distr[,i]),main = colnames(distr)[i],dotargs)
    abline(v=distr[1,i])
  }
  par(mfrow=par0$mfrow)
}

#' @export
summary.lmperm <- function(object,...){
  object$table
}


#' Method to convert into \code{Pmat} object.
#'
#'@description Convert a matrix into a \code{Pmat} object.
#'
#'@param x a matrix.
#'
#'@export
as.Pmat <- function(x) {UseMethod("as.Pmat")}

#'@export
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
as.matrix.Pmat <- function(x, ...){
  return(matrix(x, ncol = NCOL(x)))
}


#' @export
"[<-.Pmat" <- function(x,i,j,value){
  x <- as.matrix(x)
  x[i,j] <- value
  return(x)
}

#methods for np=============================
np <- function(object, ...) {UseMethod("np")}

np.matrix <-function(object){
  return(NCOL(object))
}

np.Pmat <- function(object){
  return(attr(object,which = "np"))
}

np.list <- function(object){
  return(sapply(object,function(x){np(x)}))
}






















