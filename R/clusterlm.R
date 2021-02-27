#' Cluster test for longitudinal data
#'
#' @description Compute the cluster-mass test for longitudinal linear model.
#' @param formula A formula object where the left part is a matrix defined in the global environment.
#' @param data A data frame for the independant variables.
#' @param np The number of permutations. Default value is \code{5000}.
#' @param method A character string indicating the method used to handle nuisance variables. Default is \code{NULL} and will switch to \code{"freedman_lane"} for the fixed effects model and to \code{"Rd_kheradPajouh_renaud"} for the repeated measures ANOVA. See \link{lmperm} or \link{aovperm} for details on the permutation methods.
#' @param type A character string to specify the type of transformations: "permutation" and "signflip" are available. Is overridden if P is given. See help from Pmat.
#' @param test A character string to specify the name of the test. Default is \code{"fisher"}. \code{"t"} is available for the fixed effects model.
#' @param threshold A numerical value that specify the threshold for the \code{"clustermass"} multiple comparisons procedure. If it is a vector each value will be associated to an effect. If it is scalar the same threshold will be used for each test. Default value is \code{NULL} and will compute a threshold based on the 0.95 quantile of the choosen test statistic.
#' @param aggr_FUN A function used as mass function. It should aggregate the statistics of a cluster into one scalar. Default is the sum of squares fot t statistic and sum for F statistic.
#' @param multcomp A vector of character defining the methods of multiple comparisons to compute. Default is \code{"clustermass"}, and the additional options are available : \code{"tfce"},\code{"bonferroni"}, \code{"holm"}, \code{"troendle"}, \code{"minP"} and \code{"benjamini_hochberg"}.
#' @param ... Futher arguments, see details.
#' @return A \code{clusterlm} object. Use the \link{plot.clusterlm} or \link{summary.clusterlm} method to see results of the tests.
#' @details
#' The random effects model is only avaible with a F statistic.\cr
#'
#' Other arguments could be pass in \code{...} :\cr \cr
#' \code{P} : A matrix containing the permutation of class \code{matrix} or \code{Pmat}; which is used for the reproductibility of the results. The first column must be the identity. \code{P} overwrites \code{np} argument.\cr \cr
#' \code{rnd_rotation} : A matrix of random value to compute a rotation of size \eqn{n \times n} that will be used for the \code{"huh_jhun"} method. \cr \cr
#' \code{p_scale = FALSE} : if set to \code{TRUE}, the several multiple comparisons procedures are compute on the \code{1 - p} scale, where \code{p} is the p-value. The threshold have to be set between 0 and 1 (eg: \code{threshold = 0.95}). The function \code{aggr_FUN} should be big when there is evidence against the null (eg: \code{aggr_FUN = function(p)sum(abs(log(1-p)))}. Moreover under the probability scale the cluster mass statistics is sensitive to the number permutations.\cr \cr
#' \code{H}, \code{E}, \code{ndh} : the parameters used for the \code{"tfce"} method. Default values are set to \code{H = 2} for the height parameter, to \code{E = 0.5} for the extend parameter and to \code{ndh = 500} for the number terms to approximate the integral.\cr \cr
#' \code{alpha = 0.05} : the type I error rate. Used for the \code{troendle} multiple comparisons procedure.\cr \cr
#' \code{return_distribution = FALSE} : return the permutation distribution of the statistics. Warnings : return one high dimentional matrices (number of test times number of permutation) for each test.\cr
#' \code{coding_sum} : a logical defining the coding of the design matrix to \code{contr.sum}: set by default to \code{TRUE} for ANOVA (when the argument \code{test} is \code{"fisher"} ) to tests main effects and is set to \code{FALSE} when \code{test} is \code{"t"}.  If \code{coding_sum} is set to \code{FALSE} the design matrix is computed with the coding defined in the dataframe and the tests of simple effets are possible with a coding of the dataframe set to \code{contr.treatment}. \cr
#'
#' @seealso \code{\link{plot.clusterlm}} \code{\link{summary.clusterlm}}
#'
#'@references
#'Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG-and MEG-data. Journal of neuroscience methods, 164(1), 177-190.
#'
#'Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. Neuroimage, 44(1), 83-98.
#'
#'@examples
#'
#' ## Cluster-mass for repeated measures ANOVA
#' ## Warning : np argument must be greater (recommendation: np >= 5000)
#' electrod_O1 <- clusterlm(attentionshifting_signal ~ visibility*emotion*direction
#'          + Error(id/(visibility*emotion*direction)), data = attentionshifting_design,
#'          np = 50)
#'
#' ## Results
#' plot(electrod_O1)
#'
#' ## Results with labels on the x axis that represent seconds from time-locked event:
#' plot(electrod_O1, nbbaselinepts = 200, nbptsperunit = 1024)
#'
#' ## Tables of clusters
#' electrod_O1
#'
#' \dontrun{
#' ## Change the function of the aggregation
#'
#' ## Sum of squares of F statistics
#' electrod_O1_sum <- clusterlm(attentionshifting_signal ~ visibility*emotion*direction
#'          + Error(id/(visibility*emotion*direction)), data = attentionshifting_design,
#'          aggr_FUN = function(x)sum(x^2))
#'
#' ## Length of the cluster
#' electrod_O1_length <- clusterlm(attentionshifting_signal ~ visibility*emotion*direction
#'          + Error(id/(visibility*emotion*direction)), data = attentionshifting_design,
#'          aggr_FUN = function(x)length(x))
#'
#'
#' ## All multiple comparisons procedures for repeated measures ANOVA
#' ## Permutation method "Rde_kheradPajouh_renaud"
#' full_electrod_O1 <- clusterlm(attentionshifting_signal ~ visibility*emotion*direction
#'           + Error(id/(visibility*emotion*direction)), data = attentionshifting_design,
#'           method = "Rde_kheradPajouh_renaud", multcomp = c("troendle", "tfce",
#'           "clustermass", "bonferroni", "holm", "benjamini_hochberg"))
#'}
#'
#'@author jaromil.frossard@unige.ch
#'@export
#' @family main function
clusterlm <- function(formula, data=NULL, np = 5000, method = NULL, type = "permutation", test = "fisher", threshold = NULL, aggr_FUN = NULL,
                      multcomp = "clustermass", ...){

  cl = match.call()
  if(is.null(data)){data <- model.frame(formula = formula)}



  ############
  #Formula CHECK
  Terms <- terms(formula, special = "Error", data = data)
  indError <- attr(Terms, "specials")$Error

  #dotargs
  dotargs = list(...)

  ####other parameters

  if(is.null(dotargs$alpha)){
    dotargs$alpha = 0.05
  }

  if(is.null(dotargs$p_scale)){
    dotargs$p_scale = F
  }

  if(is.null(dotargs$H)){
    switch(test,
     "t" = {dotargs$H = 2},
     "fisher" = {dotargs$H = 1})
  }

  if(is.null(dotargs$E)){
    dotargs$E = 0.5
  }


  if(is.null(dotargs$ndh)){
    dotargs$ndh = 500
  }

  #clusterdepth
  if(is.null(dotargs$border)){
    dotargs$border = "reverse"
  }
  if(is.null(dotargs$depth_scale)){
    dotargs$depth_scale = "head_and_tail"
  }


  if(is.null(dotargs$return_distribution)){
    dotargs$return_distribution = F
  }

  # if(is.null(threshold)){
  #   switch(test,
  #          "t" = {threshold = 2},
  #          "fisher" = {threshold = 4})
  # }

  if(is.null(dotargs$new_method)){
    dotargs$new_method = F
  }

  if(is.null(dotargs$coding_sum)){
    switch(test,
           "t" = {dotargs$coding_sum = F},
           "fisher" = {dotargs$coding_sum = T})
  }

  multcomp <- match.arg(multcomp, c("clustermass","clusterdepth", "tfce", "troendle","minP" , "bonferroni", "holm", "benjamini_hochberg"),
                        several.ok = T)

  ###switch fix effet
  if (is.null(indError)) {
    result <- clusterlm_fix( formula = formula, data = data, method = method, type = type, test = test, np = np,
                             P = dotargs$P, rnd_rotation = dotargs$rnd_rotation, aggr_FUN = aggr_FUN,
                             E = dotargs$E, H = dotargs$H, threshold = threshold,
                             return_distribution = dotargs$return_distribution, cl = cl, multcomp = multcomp,
                             alpha = dotargs$alpha, p_scale = dotargs$p_scale, coding_sum = dotargs$coding_sum,ndh = dotargs$ndh,
                             new_method = dotargs$new_method, border = dotargs$border, depth_scale = dotargs$depth_scale)
  } else if (!is.null(indError)){
    if(test!="fisher"){
      warning("Random effects model only accept fisher statistics. Test statistic is set to fisher.")
      test="fisher"}
    result <- clusterlm_rnd( formula = formula, data = data, method = method, type = type, test = test, np = np,
                             P = dotargs$P, rnd_rotation = dotargs$rnd_rotation, aggr_FUN = aggr_FUN,
                             E = dotargs$E, H = dotargs$H, threshold = threshold,
                             return_distribution = dotargs$return_distribution, cl = cl, multcomp = multcomp,
                             alpha = dotargs$alpha, p_scale = dotargs$p_scale, coding_sum = dotargs$coding_sum,ndh = dotargs$ndh,
                             new_method = dotargs$new_method, border = dotargs$border, depth_scale = dotargs$depth_scale)}

  ###output
  return(result)
}
