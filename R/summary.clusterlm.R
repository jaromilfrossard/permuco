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
#' @family summary
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