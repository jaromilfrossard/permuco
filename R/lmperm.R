#' Permutation tests for regression parameters
#'
#' @description Compute permutation marginal tests for linear models. This function produces t statistics with univariate and bivariate p-values. It gives the choice between multiple methods to handle nuisance variables.
#' @param formula A formula object.
#' @param data A data frame or matrix.
#' @param np The number of permutations. Default value is \code{5000}.
#' @param method A character string indicating the method use to handle nuisance variables. Default is \code{"freedman_lane"}. For the other methods, see details.
#' @param type A character string to specify the type of transformations: "permutation" and "signflip" are available. Is overridden if P is given. See help from Pmat.
#' @param ... Futher arguments, see details.
#' @return A \code{lmperm} object. See \link{aovperm}.
#' @details The following methods are available for the fixed effects model defined as \eqn{y = D\eta + X\beta + \epsilon}. If we want to test \eqn{\beta = 0} and take into account the effects of the nuisance variables \eqn{D}, we transform the data :
#' \tabular{lccc}{
#' \code{method} argument \tab \eqn{y*} \tab \eqn{D*} \tab \eqn{X*}\cr
#' \code{"draper_stoneman"} \tab \eqn{y} \tab \eqn{D} \tab \eqn{PX}\cr
#' \code{"freedman_lane"} \tab \eqn{(H_D+PR_D)y} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"manly"} \tab \eqn{Py} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"terBraak"} \tab \eqn{(H_{X,D}+PR_{X,D})y} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"kennedy"} \tab \eqn{PR_D y} \tab \tab \eqn{R_D X}\cr
#' \code{"huh_jhun"} \tab \eqn{PV'R_Dy} \tab \tab \eqn{V'R_D X}\cr
#' \code{"dekker"} \tab \eqn{y} \tab \eqn{D} \tab \eqn{PR_D X}\cr
#' }
#'
#' Other arguments could be pass in \code{...} :\cr \cr
#' \code{P} : a matrix containing the permutations of class \code{matrix} or \code{Pmat} for the reproductibility of the results. The first column must be the identity. \code{P} overwrites \code{np} argument. \cr \cr
#' \code{rnd_rotation} : a random matrix of size \eqn{n \times n} to compute the rotation used for the \code{"huh_jhun"} method.
#' @seealso \code{\link{aovperm}} \code{\link{plot.lmperm}}
#'
#' @references
#' Kherad-Pajouh, S., & Renaud, O. (2010). An exact permutation method for testing any effect in balanced and unbalanced fixed effect ANOVA. Computational Statistics & Data Analysis, 54(7), 1881-1893.
#'
#' Kherad-Pajouh, S., & Renaud, O. (2015). A general permutation approach for analyzing repeated measures ANOVA and mixed-model designs. Statistical Papers, 56(4), 947-967.
#'
#' Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M., & Nichols, T. E. (2014). Permutation inference for the general linear model. Neuroimage, 92, 381-397.
#'
#' @author jaromil.frossard@unige.ch
#' @examples
#' ## data
#' data("emergencycost")
#'
#' ## Testing at 14 days
#' emergencycost$LOS14 <- emergencycost$LOS - 14
#'
#' ## Univariate t test
#' contrasts(emergencycost$insurance) <- contr.sum
#' contrasts(emergencycost$sex) <- contr.sum
#'
#' ## Warning : np argument must be greater (recommendation: np>=5000)
#' modlm_cost_14 <- lmperm(cost ~ LOS14*sex*insurance, data = emergencycost, np = 2000)
#' modlm_cost_14
#' @export
#' @family main function
lmperm<-function(formula, data = NULL, np = 5000, method = NULL, type = "permutation",... ){
  if(is.null(data)){data <- model.frame(formula = formula)}

  ############
  #Formula CHECK
  Terms <- terms(formula, special = "Error", data = data)
  indError <- attr(Terms, "specials")$Error

  #dotargs
  dotargs=list(...)

  ###switch fix effet
  if (is.null(indError)) {
    result <- lmperm_fix( formula = formula, data = data, method = method, type = type, np = np, P = dotargs$P,
                          rnd_rotation = dotargs$rnd_rotation, new_method = dotargs$new_method)
  } else
  {
    stop("the random effects model is not implemented yet. See aovperm")
  }

  ###output
  return(result)
}
