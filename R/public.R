#' Create a set of permutations/coinflip.
#'
#' @description Compute a permutation matrix used as argument in \link{aovperm}, \link{lmperm}, \link{clusterlm} functions. The first column represents the identity permutation.
#' @param np A numeric value for the number of permutations. Default is 5000.
#' @param n A numeric value for the number of observations.
#' @param type A character string to specify the type of transformations: "permutation" and "coinflip" are available. See details.
#' @param counting A character string to specify the selection of the transformations. "all" and "random" are available. See details.
#' @return A matrix n x np containing the permutations/coinflips. First permutation is the identity.
#' @details \code{couting} can set to :\cr
#' \code{"random"} : \code{np} random with replacement permutations/coinflips among the \code{n!}/\code{2^n}  permutations.\cr
#' \code{"all"} : all \code{n!}/\code{2^n} possible permutations/coinflips.\cr
#' @importFrom permute allPerms
#' @importFrom permute how
#' @examples
#' ## data
#' data("emergencycost")
#'
#' ## Create a set of 2000 permutations
#' set.seed(42)
#' pmat = Pmat(np = 2000, n = nrow(emergencycost))
#' cfmat = Pmat(np = 2000, n = nrow(emergencycost), type = "coinflip")
#'
#' ## centrering the covariate to the mean
#' emergencycost$LOSc <- scale(emergencycost$LOS, scale = FALSE)
#'
#' ## ANCOVA
#' mod_cost_0 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost, np = 2000)
#' mod_cost_1 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost, P = pmat)
#' mod_cost_2 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost, P = pmat)
#' mod_cost_3 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost, P = cfmat)
#'
#' ## Same p-values for both models 1 and 2 but different of model 0
#' mod_cost_0
#' mod_cost_1
#' mod_cost_2
#' mod_cost_3
#'
#' @export
Pmat <- function(np = 5000, n, type = "permutation", counting = "random"){
  type <- match.arg(type, c("permutation","coinflip"))
  counting <- match.arg(counting,c("random","all"))
  #warnings type dimension
  if(type=="permutation"){
  switch(counting,
         "all" = {
           if(n > 8){
           warning("'all' permutations are not feasible for n > 8, Pmat is computed with the 'random' counting.")
             counting <- "random"}},
         {
           if(n <= 8){
             if(factorial(n) <= np){
               warning("n!<= np 'all' permutations are feasible, Pmat is computed with the 'all' counting.")
               counting <- "all"
               }
             }
           }
         )}else if(type=="coinflip")
           {switch(counting,
                  "all" = {
                    if(n > 18){
                      warning("'all' coinflip are not feasible for n > 18, Pmat is computed with the 'random' counting.")
                      counting <- "random"}},
                  {
                    if(2^n <= np){
                        warning("2^n <= np 'all' coinflip are feasible, Pmat is computed with the 'all' counting.")
                        counting <- "all"
                      }

                  }
           )}
  #creating matrix
  if(type=="permutation"){
  switch(counting,
         "random"={P <- cbind(1:n, replicate(np - 1, sample(n, n, replace = F)))},
         "all"={
             P <- t(allPerms(n,control = how(observed = T)))
             #as.matrix
             attr(P, "control") <- NULL
             attr(P, "observed") <- NULL
             class(P) <- "matrix"
             np <- factorial(n)})}
    else if(type=="coinflip"){
      switch(counting,
             "random" = {P <- cbind(rep(1, n), replicate(np - 1, sample(c(1, -1), n, replace = T)))},
             "all" = {
               P <- t(as.matrix(expand.grid(as.data.frame(t(cbind(rep(1,n),rep(-1,n)))))))
               rownames(P) <- NULL
               np <- 2^n})
             }

  attr(P,which = "type") <- type
  attr(P,which = "counting") <- counting
  attr(P,which = "np") <- np
  attr(P,which = "n") <- n
  class(P) <- "Pmat"
  return(P)
}




