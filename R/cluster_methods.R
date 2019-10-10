#' Print \code{clusterlm} object.
#'
#' @description Display the corrected p-values for each effects. Results of the \code{"clustermass"} procedure.
#'
#' @param x A \code{clusterlm} object.
#' @param multcomp A character string indicating the multiple comparison procedure to print. Default is NULL a print the first multiple comparisons procedure of the \code{clusterlm} object.
#' @param alternative A character string indicating the alternative hypothesis. Choose between \code{"two.sided"}, \code{"greater"}, \code{"less"}. Default is \code{"two.sided"}.
#' @param ... Further arguments pass to \code{print}.
#' @export
print.clusterlm <- function(x, multcomp = NULL, alternative = "two.sided", ...){
  print(summary(x, multcomp = multcomp, alternative = alternative, ... = ...))
}


#' @export
print.multcomp_table <- function(x, ...) {
  cat("Effect: ",attr(x,"effect_name"), ".\n",sep="")
  cat("Statistic: ",attr(x,"test"),"(",paste(attr(x,"df"),collapse=", "),")", ".\n",sep="")
  cat("Permutation Method: ",attr(x,"method"), ".\n",sep="")
  cat("Number of Dependant Variables: ",attr(x,"nDV"), ".\n",sep="")
  cat("Number of Permutations: ",attr(x,"np"), ".\n",sep="")
  cat("Multiple Comparisons Procedure: ",attr(x,"multcomp"), ".\n",sep="")
  if(attr(x,"multcomp") == "clustermass"){
    cat("Threshold: ",attr(x,"threshold"),".\n",sep="")
    cat("Mass Function: ",attr(x,"fun_name"),".\n",sep="")
  }
  if(attr(x,"table_type")=="cluster"){
    if(attr(x,"multcomp")=="clustermass"){
      cat("Table of clusters.\n")
    }else if(attr(x,"multcomp")!="clustermass"){
      cat("Table of pseudo-clusters.\n")
      }
    }else if(attr(x,"table_type")=="full"){
    cat("Table of all tests.\n")
  }

  cat("\n")
  if(!is.null(attr(x,"nocluster"))){
    if(attr(x,"nocluster")){
    cat("No cluster above the threshold.\n")}else{
      print.data.frame(x)}
  }else{
  print.data.frame(x)}
  cat("\n\n")
}

#' @export
print.listof_multcomp_table<- function(x,...){
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

#' Summarize of a \code{clusterlm} object.
#'
#' @description Display the corrected p-values for each effects.
#'
#' @param object A \code{clusterlm} object.
#' @param alternative A character string indicating the alternative hypothesis. Choose between \code{"two.sided"}, \code{"greater"}, \code{"less"}. Default is \code{"two.sided"}.
#' @param multcomp A character string indicating the multiple comparison procedure to display.
#' @param table_type A character string indicating the type of table to display. Choose between \code{"cluster"}, which aggregates test into pseudo-clusters (see details for the interpretations) or \code{"full"} which displays the p-values for all tests. See details for default values.
#' @param ... Further arguments see details.
#' @return A table for each effect indicating the statistics and p-values of the clusters.
#' @details It creates the full table when the number of tests is <=15 and creates a table of pseudo-clusters overwise. Note that for the \code{"troendle"} method is not based on clustering of the data and the table of pseudo-clusters should only be used to facilitate the reading of the results.
#' @export
summary.clusterlm <- function(object, alternative = "two.sided", multcomp = NULL, table_type = NULL, ...){
  dotargs = list(...)
  if(is.null(table_type)){
    if(ncol(object$y)>15){
      table_type = "cluster"
    }else{
      table_type = "full"
    }
  }

  if(is.null(multcomp)){multcomp <- object$multcomp[1]}

  switch(table_type,
    "cluster" = {getTable = function(x, multcomp,... ){cluster_table(x, multcomp = multcomp, ... = ...)}},
    "full" = {getTable = function(x, multcomp,... ){full_table(x, multcomp = multcomp, ... = ...)}}
  )

  switch(alternative,
           "two.sided" = {
             return(getTable(object$multiple_comparison, multcomp = multcomp, ... = ...))},
           "greater" = {
             return(getTable(object$multiple_comparison_greater, multcomp = multcomp, ... = ...))},
           "less" = {
             return(getTable(object$multiple_comparison_less, multcomp = multcomp, ... = ...))})

}

# #' @export
# summary.cluster_table <- function(object,...){
#   object
# }


# summary_multcomp <- function(object, multcomp, alternative){
#   if(!(multcomp %in% object$multcomp)){
#     stop(paste("The choosen multiple comparisons procedure does not match with the ones computed in the object.
#                   Choose between ", paste(object$multcomp,sep=" "),".", sep=""))
#   }
#   switch(alternative,
#          "two.sided" = {mc = object$multiple_comparison},
#          "greater" = {mc = object$multiple_comparison_greater},
#          "less" = {mc = object$multiple_comparison_less})
#
#   out = lapply(seq_along(mc),function(i){
#     out = mc[[i]][[multcomp]]$main[,1:2,drop = F]
#     colnames(out) = paste(names(mc)[i], colnames(out))
#     out})
#   out = do.call("cbind",out)
#   return(out)
# }



#'Plot cluster or parameters.
#'
#' @description Plot method for class \code{clusterlm}.
#'
#' @param x A \code{clusterlm} object.
#' @param effect A vector of character naming the effects to display. Default is \code{"all"}.
#' @param type A character string that specified the values to highlight. \code{"statistic"} or \code{"coef"} are available. Default is \code{"statistic"}.
#' @param multcomp A character string specifying the method use to correct the p-value. It should match the one computed in the object. Default is the (first) method in the call to \link{clusterlm}. See \link{clusterlm}.
#' @param alternative A character string specifying the alternative hypothesis for the t-test. The available options are \code{"greater"}, \code{"less"} and \code{"two.sided"}. Default is \code{"two.sided"}.
#' @param enhanced_stat A logical. Default is \code{FALSE}. If \code{TRUE}, the enhanced statistic will be plotted overwise it will plot the observed statistic. Change for the \code{"tfce"} or the \code{"clustermass"} multiple comparisons procedures.
#' @param nbbaselinepts An integer. Default is 0. If the origin of the x axis should be shifted to show the start of the time lock, provide the number of baseline time points.
#' @param nbptsperunit An integer. Default is 1. Modify this value to change the scale of the label from the number of points to the desired unit. If points are e.g. sampled at 1024Hz, set to 1024 to scale into seconds and to 1.024 to scale into milliseconds.
#' @param distinctDVs Boolean. Should the DVs be plotted distictively, i.e. should the points be unlinked and should the name of the DVs be printed on the x axis ? Default is FALSE if the number of DV is large thant 15 or if the method is "clustermass" or "tfce".
#' @param ... further argument pass to plot.
#' @importFrom graphics points axis
#' @export
plot.clusterlm <- function(x, effect = "all", type = "statistic", multcomp = x$multcomp[1], alternative = "two.sided", enhanced_stat = FALSE,
                           nbbaselinepts=0, nbptsperunit=1, distinctDVs=NULL, ...) {

  ##select effect
  if("all" %in% effect){effect = names(x$multiple_comparison)}
  else if(sum(names(x$multiple_comparison)%in%effect) == 0){
    warning(" the specified effects do not exist. Plot 'all' effects.")
    effect = names(x$multiple_comparison)
  }
  effect_sel <- names(x$multiple_comparison)%in%effect

  ###switch mult comp
  switch(alternative,
         "two.sided" = {multiple_comparison = x$multiple_comparison[effect_sel]},
         "greater" = {multiple_comparison = x$multiple_comparison_greater[effect_sel]},
         "less" = {multiple_comparison = x$multiple_comparison_less[effect_sel]})

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
                      switch (alternative,
                              "less" ={hl <- -c(abs(x$threshold))},
                              "greater" ={hl <- c(abs(x$threshold))},
                              "two.sided" ={hl <- c(abs(x$threshold))}
                      )})}
         })

  title =paste(title," : ", multcomp, " correction",sep="", collapse = "")



  #####
  p = sum(NROW(data))
  rnames = row.names(data)
  cnames = colnames(data)
  nbDV = ncol(data)



  ### should the DVs be named in plot and separated ?
  if (is.null(distinctDVs)){
    if (multcomp %in% c("clustermass", "tfce"))
      distinctDVs = FALSE
    else distinctDVs = (nbDV<16)
  }
  if ((distinctDVs==TRUE) &&  (multcomp %in% c("clustermass", "tfce")))
    warning("Computations and corrections have been based on adjacency of DVs but the the plot will show separated DVs")
  #####

  par0 <- list(mfcol = par()$mfcol,mar = par()$mar,oma = par()$oma)

  par(mfcol = c(p,1),mar = c(0,4,0,0),oma = c(4,0,4,1),...=...)
  for (i in 1:p) {
    if (distinctDVs) {
      plot((1:ncol(data)-nbbaselinepts)/nbptsperunit,
           data[i,],type = "p", xaxt = "n",xlab = "",ylab = rnames[i], pch=18, cex=2,
      )
      if(i==p) axis(1, at= (1:ncol(data)-nbbaselinepts)/nbptsperunit, label=cnames)
    }
    else{
      if(i==p){xaxt = NULL}else{xaxt = "n"}
      plot((1:ncol(data)-nbbaselinepts)/nbptsperunit,
           data[i,],type = "l", xaxt = xaxt,xlab = "",ylab = rnames[i]
      )
    }
    if(type == "statistic"){
      xi = which(pvalue[i,]< x$alpha)
      y = data[i,xi]
      col="red"
      #lines(x = x,y= y,lwd=par()$lwd*2,col=col)
      points(x = (xi-nbbaselinepts)/nbptsperunit, y = y, pch=18,col=col, cex=distinctDVs+1)
      if(multcomp=="clustermass"){
        abline(h=hl[i],lty=3)
        if(x$test=="t"&alternative=="two.sided"){
          abline(h=-hl[i],lty=3)
        }
      }
    }}
  title(title,outer = T,cex = 2)
  par(mfcol = par0$mfcol, mar = par0$mar, oma = par0$oma)
}



