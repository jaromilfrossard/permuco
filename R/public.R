#' Create a set of permutations.
#'
#' @description Compute a permutation matrix used as argument in \link{aovperm}, \link{lmperm}, \link{clusterlm} functions. The first column represents the identity permutation.
#' @param np A numeric value for the number of permutations. Default is 5000.
#' @param n A numeric value for the number of observations.
#' @param type A character string to specify the type of matrix. See Details.
#' @return A matrix n x np containing the permutations/coinflips. First permutation is the identity.
#' @details \code{type} can set to :\cr
#' \code{"default"} : \code{np} random with replacement permutations among the \code{n!} permutations.\cr
# \code{"unique"} : \code{np} random without replacement permutations among the \code{n!} permutations. Not available.\cr
#' \code{"all"} : all \code{n!} possible permutations.\cr
# \code{"coinflip"} : \code{np} coinflips vectors.
#' @importFrom permute allPerms
#' @importFrom permute how
#' @examples \dontrun{
#'
#' ## Create a set of 5000 the permutation
#' set.seed(42)
#' pmat = Pmat(np = 5000, n = nrow(emergencycost))
#'
#' ## centrering the covariate to the mean
#' emergencycost$LOSc <- scale(emergencycost$LOS, scale = F)
#'
#' ## ANCOVA
#' mod_cost_0 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost)
#' mod_cost_1 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost, P = pmat)
#' mod_cost_2 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost, P = pmat)
#'
#' ## Identical p-values for both models 1 and 2 but differents of model 0
#' mod_cost_0
#' mod_cost_1
#' mod_cost_2
#'
#' }
#' @export
Pmat <- function(np = 5000, n, type = "default"){
  #warnings
  switch(type,
         "all" = {
           if(n > 8){
           warning("all type is not feasible for n > 8, Pmat is computed with the default type.")
           type <- "default"}},
         {
           if(n <= 8){
             if(factorial(n) <= np){
               warning("n!<= np all permutations are feasible, Pmat is computed with the 'all' type.")
               type <- "all"
               }
             }
           }
         )
  #matrix
  switch(type,
         "default"={P <- cbind(1:n, replicate(np - 1, sample(n, n, replace = F)))},
         "unique"= {
           P <- cbind(1:n, replicate(np - 1, sample(n, n, replace = F)))
           type <- "default"
           warning("unique type is not implemented, Pmat is computed with the 'default' type.")
         },
         "all"={
             P <- t(allPerms(n,control = how(observed = T)))
             #as.matrix
             attr(P, "control") <- NULL
             attr(P, "observed") <- NULL
             class(P) <- "matrix"
             np <- factorial(n)},
         "coinflip" = {P <- cbind(rep(1, n), replicate(np - 1, sample(c(1, -1), n, replace = T)))})
  attr(P,which = "type") <- type
  attr(P,which = "np") <- np
  attr(P,which = "n") <- n
  class(P) <- "Pmat"
  return(P)
}




