############################################################################################################
#' Cluster-depth correction (from the head only)
#'
#' @description Compute the cluster-depth test correction (from the head) given a matrix a permuted statistical signals.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param threshold A scalar that represents the threshold to create the clusters.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @export
#' @family multcomp
compute_clusterdepth_head <- function(distribution, threshold, alternative = "two.sided"){


  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

  cluster <- get_cluster(distribution = distribution, threshold = threshold, alternative = alternative, side = "starting")
  depth_head <- get_clusterdepth_head(cluster, border = "ignore")
  distr_head <- depth_distribution(distribution, head_mat = depth_head)



  pvalue <- rep(NA,ncol(cluster))
  max_cl_size <- max(table(cluster[1,cluster[1,]!=0]))


  for(cli  in seq_len(max(cluster[1,]))){
    sample <- which(cluster[1,]==cli)
    di <- distr_head#[,seq_len(max_cl_size), drop=FALSE]
    di <- rbind(c(distribution[1, sample],rep(0,ncol(di)-length(sample))),di)
    pvalue[sample] <- compute_troendle(di,
                                       alternative = alternative)$main[seq_along(sample),2]


  }

  list(main =cbind(statistic = distribution[1,], pvalue = pvalue,cluster_id = cluster[1,]),threshold = threshold)


}
#
#
#
# compute_clusterdepth_tail <- function(distribution, threshold, type ="tail", border = "keep"){
#   type <- match.arg(type,choices = c("tail", "depth"))
#
#   border <- match.arg(border,choices = c("keep","rev","rm"))
#   cluster <- compute_cluster_cpp(distribution = distribution, threshold = threshold)
#
#
#   if(type %in% c("tail")){
#     depth_tail <- cluster_depth_tail(cluster, border = border)
#     if(border=="rm"){cluster <- compute_cluster_cpp(depth_tail,0.5)}
#
#     distr_tail <- distribution_depth(distribution, tail_mat = depth_tail)}
#
#   if(type %in% c("depth")){
#     depth_head <- cluster_depth_head(cluster,border = border)
#     depth_tail <- cluster_depth_tail(cluster,border = border)
#     distr_tail <- distribution_depth(distribution, head_mat = depth_head, tail_mat = depth_tail)
#     distr_tail <- distr_tail[,rev(seq_len(ncol(distr_tail)))]}
#
#
#
#   pvalue = rep(NA,ncol(cluster))
#
#
#   for(cli  in seq_len(max(cluster[1,]))){
#     sample <- which(cluster[1,]==cli)
#     stats <- distribution[1,sample]
#     pvalue[sample] <- compute_troendle(rbind(stats,distr_tail[,seq_along(sample)+ncol(distr_tail)-length(sample),drop=F]),
#                                        alternative = "greater")$main[,2]
#
#
#   }
#
#   list(main =cbind(statistic = distribution[1,], pvalue = pvalue))
#
#
# }


