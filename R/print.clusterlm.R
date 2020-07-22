#' Print \code{clusterlm} object.
#'
#' @description Display the corrected p-values for each effects. Results of the \code{"clustermass"} procedure.
#'
#' @param x A \code{clusterlm} object.
#' @param multcomp A character string indicating the multiple comparison procedure to print. Default is NULL a print the first multiple comparisons procedure of the \code{clusterlm} object.
#' @param alternative A character string indicating the alternative hypothesis. Choose between \code{"two.sided"}, \code{"greater"}, \code{"less"}. Default is \code{"two.sided"}.
#' @param ... Further arguments pass to \code{print}.
#' @export
#' @family summary
print.clusterlm <- function(x, multcomp = NULL, alternative = "two.sided", ...){
  print(summary(x, multcomp = multcomp, alternative = alternative, ... = ...))
}