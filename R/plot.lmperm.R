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
#'@family plot
plot.lmperm <- function(x, FUN = density, ...){

  par0 <-  par()
  #data
  distr <- x$distribution


  dotargs <- list(...)

  if(is.null(dotargs$effect)){
    effect <- colnames(distr)
  }else{effect <- colnames(distr)[which(colnames(distr)%in%dotargs$effect)]}

  dotargs_par <- dotargs[names(dotargs)%in%names(par())]
  dotargs <-  dotargs[!names(dotargs)%in%c(names(par()),"effect")]




  distr <- distr[,which(colnames(distr)%in%effect),drop=F]

  #subplot
  p <- NCOL(distr)
  div <- seq_len(abs(p))
  factors <- div[p %% div == 0L]
  mfrow1 <- factors[ceiling(length(factors)/2)]
  mfrow <- c(mfrow1,p/mfrow1)

  ### param default
  if(is.null(dotargs_par$mfrow)){dotargs_par$mfrow = mfrow}
  par(dotargs_par)


  #plot
  for(i in 1:NCOL(distr)){
    argi = dotargs
    argi$x = FUN(distr[,i])
    argi$main = colnames(distr)[i]
    do.call("plot",argi)
    #plot(FUN(distr[,i]),main = colnames(distr)[i])
    abline(v=distr[1,i])
  }
  par0 <- par0[!names(par0)%in%c("cin","cra","csi","cxy","din","page")]
  par(par0)
}