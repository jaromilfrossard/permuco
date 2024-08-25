
#methods for np=============================
np <- function(object) {UseMethod("np")}

np.matrix <-function(object){
  return(NCOL(object))
}

np.Pmat <- function(object){
  return(attr(object,which = "np"))
}

np.list <- function(object){
  return(sapply(object,function(x){np(x)}))
}





