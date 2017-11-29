#' Print \code{clusterlm} object.
#'
#' @description Display with the corrected p-values for each effects. Results of the \code{"clustermass"} procedure.
#'
#' @param x A \code{clusterlm} object.
#' @param laterality A character string indicating the laterality of the tests. Choose between \code{"bilateral"}, \code{"right"}, \code{"left"}. Default is \code{"bilateral"}.
#' @param ... Further arguments pass to \code{print}.
#' @export
print.clusterlm <- function(x, laterality = "bilateral",...){
  cat(
    "Cluster ", x$test,
    " test using ", x$method,
    " to handle nuisance variables \n with ", paste(np(x$P), sep= ", ", collapse = ", "),
    " permutations and ",x$fun_name," as mass function.\n\n", sep = "")
  cat("Alternative Hypothesis : ",laterality,".\n \n",sep = "")
  switch(laterality,
         "bilateral" = {
           print(x$cluster_table,...)},
         "right" = {
           print(x$cluster_table_right,...)},
         "left" = {
           print(x$cluster_table_left,...)})
  }


#' @export
print.cluster_table <- function(x, ...) {
  cat(attr(x,"effect_name"),", threshold = ",attr(x,"threshold"),".\n",sep="")
  print.data.frame(x)
  cat("\n")
}

#' @export
print.listof_cluster_table<- function(x,...){
  dotargs = list(...)
  ei = NULL
  if(!is.null(dotargs$effect)){
    ei = which(names(x)%in%dotargs$effect)}
  if(length(ei)>0){
    for(i in ei){print(x[[i]])}
    }else{
      for(i in 1:length(x)){
        print(x[[i]])}}
}

###########################summary

#' Summarize a \code{clusterlm} object.
#'
#' @description Display the clusters with the corrected p-values for each effects. Results of the \code{"clustermass"} procedure.
#'
#' @param object A \code{clusterlm} object.
#' @param laterality A character string indicating the laterality of the tests. Choose between \code{"bilateral"}, \code{"right"}, \code{"left"}. Default is \code{"bilateral"}.
#' @param ... Further arguments see details.
#' @return A table for each effect indicating the statistics and p-values of the clusters.
#' @details If the \code{multcomp} argument is a character string that matches the \code{multcomp} argument of the \code{clusterlm} object, this method returns a matrix with the corrected statistics and p-values in columns and multiple tests by rows.
#' @export
summary.clusterlm <- function(object, laterality = "bilateral",...){
  dotargs = list(...)
  if(is.null(dotargs$multcomp)){
    switch(laterality,
           "bilateral" = {
             return(object$cluster_table)},
           "right" = {
             return(object$cluster_table_right)},
           "left" = {
             return(object$cluster_table_left)})
    }else{summary_multcomp(object = object, multcomp = dotargs$multcomp, laterality = laterality)}
}

#' @export
summary.cluster_table <- function(object,...){
  object
}


summary_multcomp <- function(object, multcomp, laterality){
  if(!(multcomp %in% object$multcomp)){
    stop(paste("The choosen multiple comparisons procedure does not match with the ones computed in the object.
                  Choose between ", paste(object$multcomp,sep=" "),".", sep=""))
  }
  switch(laterality,
         "bilateral" = {mc = object$multiple_comparison},
         "right" = {mc = object$multiple_comparison_right},
         "left" = {mc = object$multiple_comparison_left})

  out = lapply(seq_along(mc),function(i){
    out = mc[[i]][[multcomp]]$main[,1:2,drop = F]
    colnames(out) = paste(names(mc)[i], colnames(out))
    out})
  out = do.call("cbind",out)
  return(out)
}



#'Plot cluster or parameters.
#'
#' @description Plot method for class \code{clusterlm}.
#'
#' @param x A \code{clusterlm} object.
#' @param effect A vector of character naming the effects to display. Default is \code{"all"}.
#' @param type A character string that specified the values to highlight. \code{"statistic"} or \code{"coef"} are available. Default is \code{"statistic"}.
#' @param multcomp A character sting specifying the p-value to plot. Default is \code{"clustermass"}. See \link{clusterlm}.
#' @param laterality A character string specifying the laterality of the test when t-test are computed. Avaible options are \code{"right"}, \code{"left"} and \code{"bilateral"}. Default is \code{"bilateral"}.
#' @param enhanced_stat logical. Default is \code{F}. If \code{TRUE}, the enhanced statistic will be plotted overwise it will plot the observed statistic. Change for the \code{"tfce"} or the \code{"clustermass"}, multiple comparisons method.
#' @param ... further argument pass to plot.
#' @importFrom graphics points axis
#' @export
plot.clusterlm <- function(x, effect = "all", type = "statistic", multcomp = "clustermass", laterality = "bilateral", enhanced_stat = F,...) {

  ##select effect
  if("all" %in% effect){effect = names(x$multiple_comparison)}
  else if((names(x$multiple_comparison)%in%effect) == 0){
    warning(" the specified effects do not exist. Plot 'all' effects.")
    effect = names(x$multiple_comparison)
  }
  effect_sel <- names(x$multiple_comparison)%in%effect

  ###switch mult comp
switch(laterality,
       "bilateral" = {multiple_comparison = x$multiple_comparison[effect_sel]},
       "right" = {multiple_comparison = x$multiple_comparison_right[effect_sel]},
       "left" = {multiple_comparison = x$multiple_comparison_left[effect_sel]})

  pvalue = t(sapply(multiple_comparison,function(m){
    m[[multcomp]]$main[,2]}))

  statistic = t(sapply(multiple_comparison,function(m){
    m[["uncorrected"]]$main[,1]}))
  if(enhanced_stat){
    statistic = t(sapply(multiple_comparison,function(m){
      m[[multcomp]]$main[,1]}))
  }



  ##swich value
  switch(type,
         "coef"={
           data <- x$coef[effect_sel,]
           title <- "coefficients"
           hl <- NULL
         },
         "statistic" ={
           data <- statistic
           title <- paste(x$test, " statistic",sep="",collapse = "")
           if(multcomp=="clustermass"){
           switch(x$test,
                  "fisher"={hl <- x$threshold},
                  "t"={
                    switch (laterality,
                      "left" ={hl <- c(-abs(x$threshold))},
                      "right" ={hl <- c(abs(x$threshold))},
                      "bilateral" ={hl <- c(-x$threshold,x$threshold)}
                    )})}
         })

  title =paste(title," : ", multcomp, " p-values",sep="", collapse = "")


  #####
  p = sum(NROW(data))
  rnames = row.names(data)


  #####

  par0 <- list(mfcol = par()$mfcol,mar = par()$mar,oma = par()$oma)

  par(mfcol = c(p,1),mar = c(0,4,0,0),oma = c(4,0,4,1),...=...)
  for (i in 1:p) {
    if(i==p){xaxt = NULL}else{xaxt = "n"}
    plot(
      data[i,],type = "l",xaxt = xaxt,xlab = "",ylab = rnames[i]
    )
    if(type == "statistic"){
      xi = which(pvalue[i,]< x$alpha)
      y = data[i,xi]
      col="red"
      #lines(x = x,y= y,lwd=par()$lwd*2,col=col)
      points(x = xi, y = y, pch=par()$pch,col=col)
      if(multcomp=="clustermass"){abline(h=hl[i],lty=3)}
  }}
  title(title,outer = T,cex = 2)
  par(mfcol = par0$mfcol, mar = par0$mar, oma = par0$oma)
  }


