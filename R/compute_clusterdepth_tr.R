#' ############################################################################################################
#' #' Cluster-depth correction (head-tail combined)
#' #'
#' #' @description Compute the clusterdepth test correction given a matrix a permuted statistical signals.
#' #'
#' #' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' #' @param threshold A scalar that represents the threshold to create the clusters.
#' #' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' #' @details The final troendle correction is computed by taking all tests (head and tail) together.
#' #' @export
#' #' @family multcomp
#' compute_clusterdepth_tr <- function(distribution,threshold, alternative = "two.sided"){
#'   alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
#'
#'
#'   cluster <- get_cluster(distribution = distribution[1,,drop=F],
#'                          threshold = threshold, alternative = alternative,side= "all")
#'
#'   cluster_head <- get_cluster(distribution = distribution, threshold = threshold, alternative = alternative, side = "starting")
#'   depth_head <- get_clusterdepth_head(cluster_head, border = "ignore")
#'   distr_head <- depth_distribution(distribution, head_mat = depth_head)
#'
#'
#'   cluster_tail <- get_cluster(distribution = distribution, threshold = threshold, alternative = alternative, side = "ending")
#'   depth_tail <- get_clusterdepth_tail(cluster_tail,border = "ignore")
#'   distr_tail <- depth_distribution(distribution, tail_mat = depth_tail)
#'
#'
#'   pvalues <- rep(NA,ncol(cluster_tail))
#'
#'
#'
#'   cl_id <- unique(cluster[(cluster_head[1,]!=0)&(cluster_tail[1,]!=0)])
#'
#'
#'   for(cli  in cl_id){
#'     sample <- which(cluster[1,]==cli)
#'     dih <- distr_head
#'     dih <- rbind(c(distribution[1, sample],rep(0,ncol(dih)-length(sample))),dih)
#'
#'
#'     dit <- distr_tail
#'     dit <- rbind(c(rep(0,ncol(dit)-length(sample)), distribution[1, sample]),
#'                 dit)
#'     pvi <- compute_troendle(cbind(dih,dit), alternative = alternative)$main[,2]
#'
#'     pvalues[sample] = pmax(pvi[seq_along(sample)],
#'                             pvi[seq_along(sample)+ncol(dih)+ ncol(dit)-length(sample)])
#'
#'     # stats <- distribution[1,sample]
#'     # stats <- c(rep(0,max_cl_size-length(stats)),stats)
#'     # pvalue_tail[sample] <- compute_troendle(rbind(stats,distr_tail[,seq_len(max_cl_size),drop=F]),
#'     #                                         alternative = alternative)$main[seq_along(sample),2]
#'
#'
#'   }
#'
#'
#'   list(main =cbind(statistic = distribution[1,], pvalue = pvalues, cluster_id = as.integer(cluster)),threshold = threshold)
#'
#'
#' }
#'
#'
#'
