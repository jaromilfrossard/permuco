#' P-values based on permutation tests for ANOVA and repeated measures ANOVA designs.
#'
#' @description  Provides p-values for omnibus tests based on permutations for factorial and repeated measures ANOVA. This function produces the F statistics, the parametric, p-values based on Gaussian and sphericity assumptions and p-values based on the permutation method that handle nuisance variables.
#' @param formula A formula object. The formula for repeated measures ANOVA should be written using the same notation as \link{aov} by adding
#' \code{+Error(id/within)}, where \code{id} is the factor that identify the subjects and \code{within} is the within factors.
#' @param data A data frame or matrix.
#' @param np The number of permutations. Default value is \code{5000}.
#' @param method A character string indicating the method used to handle nuisance variables. Default is \code{NULL} and will change if to \code{"freedman_lane"} for the fixed effects model and \code{"Rd_kheradPajouh_renaud"} for the random effect models. See Details for other methods.
#' @param ... Futher arguments, see details.
#'
#' @return A \code{lmperm} object containing most of the objects given in an \link{lm} object, an ANOVA table with parametric and permutation p-values, the test statistics and the permutation distributions.
#'
#' @details The following methods are available for the fixed effects model defined as \eqn{y = D\eta + X\beta + \epsilon}. If we want to test \eqn{\beta = 0} and take into account the effects of the nuisance variables \eqn{D}, we transform the data :
#' \tabular{lccc}{
#' \code{method} argument \tab \eqn{y*} \tab \eqn{D*} \tab \eqn{X*}\cr
#' \code{"draper_stoneman"} \tab \eqn{y} \tab \eqn{D} \tab \eqn{PX}\cr
#' \code{"freedman_lane"} \tab \eqn{(H_D+PR_D)y} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"manly"} \tab \eqn{Py} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"terBraak"} \tab \eqn{(H_{X,D}+PR_{X,D})y} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"kennedy"} \tab \eqn{PR_D y} \tab \tab \eqn{R_D X}\cr
#' \code{"huh_jhun"} \tab \eqn{PV'R_Dy} \tab \tab \eqn{PV'R_D X}\cr
#' \code{"dekker"} \tab \eqn{y} \tab \eqn{D} \tab \eqn{PR_D X}\cr
#' }
#'The following methods are available for the random effects model \eqn{y = D\eta + X\beta + E\kappa + Z\gamma+ \epsilon}. If we want to test \eqn{\beta = 0} and take into account the effect of the nuisance variable \eqn{D} we can transform the data by permutation:
#' \tabular{lccccc}{
#' \code{method} argument \tab \eqn{y*} \tab \eqn{D*} \tab \eqn{X*} \tab \eqn{E*} \tab \eqn{Z*}\cr
#' \code{"Rd_kheradPajouh_renaud"} \tab \eqn{R_D y} \tab  \tab \eqn{R_D X} \tab \eqn{R_D E} \tab \eqn{R_D Z}\cr
#' \code{"Rde_kheradPajouh_renaud"} \tab \eqn{R_{D,E}y}  \tab  \tab \eqn{R_{D,E} X} \tab   \tab \eqn{R_{D,E}Z}\cr}
#'
#' Other arguments could be pass in \code{...} :\cr \cr
#' \code{P} : a matrix,  of class \code{matrix} or \code{Pmat}, containing the permutations (for the reproductibility of the results). The first column must be the identity permutation (not checked). \code{P} overwrites \code{np} argument. \cr
#' \code{rnd_rotation} : a random matrix of size \eqn{n \times n} to compute the rotation used for the \code{"huh_jhun"} method.
#' \code{coding_sum} : a logical set to \code{TRUE} defining the coding of the design matrix to \code{contr.sum} to test the main effects. If it is set to \code{FALSE} the design matrix is computed with the coding defined in the dataframe. The tests of simple effets are possible with a coding of the factors of the dataframe set to \code{contr.treatment}.\cr
#' @seealso \code{\link{lmperm}} \code{\link{plot.lmperm}}
#' @author jaromil.frossard@unige.ch
#' @importFrom stats terms contr.sum
#'
#' @examples
#' ## data
#' data("emergencycost")
#'
#' ## centrering the covariate to the mean
#' emergencycost$LOSc <- scale(emergencycost$LOS, scale = FALSE)
#'
#' ## ANCOVA
#' mod_cost_0 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost, np = 2000)
#' mod_cost_0
#'
#' ## Testing at 14 days
#' emergencycost$LOS14 <- emergencycost$LOS - 14
#'
#' mod_cost_14 <- aovperm(cost ~ LOS14*sex*insurance, data = emergencycost, np = 2000)
#' mod_cost_14
#'
#' ## Effect of sex within the public insured
#' contrasts(emergencycost$insurance) <- contr.treatment
#' contrasts(emergencycost$sex) <- contr.sum
#' emergencycost$insurance <- relevel(emergencycost$insurance, ref = "public")
#'
#' mod_cost_se <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost,
#'                         np = 2000, coding_sum = FALSE)
#' mod_cost_se
#'
#'
#' ## Repeated measures ANCOVA
#' ## data
#' data(jpah2016)
#'
#' ## centrering the covariate
#' jpah2016$bmic <- scale(jpah2016$bmi, scale = FALSE)
#'
#' mod_jpah2016 <- aovperm(iapa ~ bmic*condition*time+ Error(id/(time)),
#'                     data = jpah2016, method = "Rd_kheradPajouh_renaud")
#' mod_jpah2016
#'
#'
#' @export
aovperm<-function(formula, data=NULL, np = 5000, method = NULL,...){
  #method <- pmatch(method)

  if(is.null(data)){data <- model.frame(formula = formula)}

  #Formula CHECK
  Terms <- terms(formula, special = "Error", data = data)
  indError <- attr(Terms, "specials")$Error

  #check for intercept
  if(!attr(Terms,"intercept")){warning("Intercept should be specified in the formula")}

  #dotargs
  dotargs=list(...)

  ###switch fix effet
  if (is.null(indError)) {
    result <- aovperm_fix( formula = formula, data = data, method = method, np = np, coding_sum = dotargs$coding_sum, P = dotargs$P,
                           rnd_rotation = dotargs$rnd_rotation, new_method = dotargs$new_method)
  } else if (!is.null(indError))
  {
    result <- aovperm_rnd( formula = formula, data = data, method = method, np = np, coding_sum = dotargs$coding_sum, P = dotargs$P,
                           rnd_rotation = dotargs$rnd_rotation, new_method = dotargs$new_method)
  }

  ###output
  return(result)
}

#' Permutation test for regression parameters
#'
#' @description Compute permutation marginal test for linear model. This function produces t statistics with univariate and bivariate p-values. It gives the choice between multiple methods to handle nuisance variables.
#' @param formula A formula object.
#' @param data A data frame or matrix.
#' @param np The number of permutations. Default value is \code{5000}.
#' @param method A character string indicating the method use to handle nuisance variables. Default is \code{"freedman_lane"}. Se details for the other methods.
#' @param ... Futher arguments, see details.
#' @return A \code{lmperm} object. see \link{aovperm}.
#' @details The following methods are available for the fixed effects model defined as \eqn{y = D\eta + X\beta + \epsilon}. If we want to test \eqn{\beta = 0} and take into account the effects of the nuisance variables \eqn{D}, we transform the data :
#' \tabular{lccc}{
#' \code{method} argument \tab \eqn{y} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"draper_stoneman"} \tab \eqn{y} \tab \eqn{D} \tab \eqn{PX}\cr
#' \code{"freedman_lane"} \tab \eqn{(H_D+PR_D)y} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"manly"} \tab \eqn{Py} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"terBraak"} \tab \eqn{(H_{X,D}+PR_{X,D})y} \tab \eqn{D} \tab \eqn{X}\cr
#' \code{"kennedy"} \tab \eqn{PR_D y} \tab \tab \eqn{R_D X}\cr
#' \code{"huh_jhun"} \tab \eqn{PV'R_Dy} \tab \tab \eqn{PV'R_D X}\cr
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
#' modlm_cost_14 <- lmperm(cost ~ LOS14*sex*insurance, data = emergencycost, np = 2000)
#' modlm_cost_14
#' @export
lmperm<-function(formula, data = NULL, np = 5000, method = NULL,... ){
  if(is.null(data)){data <- model.frame(formula = formula)}

  ############
  #Formula CHECK
  Terms <- terms(formula, special = "Error", data = data)
  indError <- attr(Terms, "specials")$Error

  #dotargs
  dotargs=list(...)

  ###switch fix effet
  if (is.null(indError)) {
    result <- lmperm_fix( formula = formula, data = data, method = method, np = np, P = dotargs$P,
                          rnd_rotation = dotargs$rnd_rotation, new_method = dotargs$new_method)
  } else
  {
    stop("the random effects model is not implemented yet.")
  }

  ###output
  return(result)
}
