############################################################################################################
#' Clusterdept correction
#'
#' @description Compute the clusterdepth test correction given a matrix a permuted statistical signals.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param threshold A scalar that represents the threshold to create the clusters.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @param depth_scale A character string indicating the measure of the depth. Default is \code{"head_and_tail"} corrects independently from head and tail. \code{"unique_depth"} is also available and use a unique depth measure.
#' @param border A character string indicating the method to handle the clusters at the border. Default is \code{"reverse"} and will reverse the measure of the depth. \code{"ignore"} is also available and consider these cluster similarly to the others.
#' @details The \code{aggr_FUN} argument may take predefined function as the sum: \code{aggr_FUN = sum} and also user-defined function as the sum of squares: \code{aggr_FUN = function(x){sum(x^2)}}
#' @export
#' @family multcomp
compute_clusterdepth <- function(distribution,threshold, alternative = "two.sided", depth_scale = "head_and_tail", border = "reverse"){
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  depth_scale <- match.arg(depth_scale , choices = c("unique_depth", "head_and_tail"))
  border <- match.arg(border, choices = c("reverse", "ignore"))


  cluster <- get_cluster(distribution = distribution, threshold = threshold, alternative = alternative)
  depth_head <- get_clusterdepth_head(cluster, border = border)
  depth_tail <- get_clusterdepth_tail(cluster,border = border)

  if(depth_scale %in% c("head_and_tail")){
    distr_head <- depth_distribution(distribution, head_mat = depth_head)
    distr_tail <- depth_distribution(distribution, tail_mat = depth_tail)
  }

  if(depth_scale %in% c("unique_depth")){
    distr_head <- depth_distribution(distribution, head_mat = depth_head, tail_mat = depth_tail)
    distr_tail <- distr_head[,rev(seq_len(ncol(distr_head)))]
  }





  pvalue = rep(NA,ncol(cluster))


  for(cli  in seq_len(max(cluster[1,]))){
    sample <- which(cluster[1,]==cli)
    stats <- distribution[1,sample]
    headp <- compute_troendle(rbind(stats,distr_head[,seq_along(sample),drop=F]),
                              alternative = alternative)$main[,2]
    tailp <- compute_troendle(rbind(stats,distr_tail[,seq_along(sample)+ncol(distr_tail)-length(sample),drop=F]),
                              alternative = alternative)$main[,2]
    pvalue[sample] <- apply(cbind(headp,tailp),1,max)


  }

  list(main =cbind(statistic = distribution[1,], pvalue = pvalue),depth_scale = depth_scale, border = border)


}





# compute_clusterdepth_head <- function(distribution, threshold, type ="head", border = "keep"){
#   type <- match.arg(type,choices = c("head", "depth"))
#
#   border <- match.arg(border,choices = c("keep","rev","rm"))
#   cluster <- compute_cluster_cpp(distribution = distribution, threshold = threshold)
#
#
#   if(type %in% c("head")){
#     depth_head <- cluster_depth_head(cluster, border = border)
#     distr_head <- distribution_depth(distribution, head_mat = depth_head)
#
#     if(border=="rm"){ cluster <- compute_cluster_cpp( depth_head,0.5)}
#   }
#
#   if(type %in% c("depth")){
#     depth_head <- cluster_depth_head(cluster,border = border)
#     depth_tail <- cluster_depth_tail(cluster,border = border)
#     distr_head <- distribution_depth(distribution, head_mat = depth_head, tail_mat = depth_tail)}
#
#
#
#   pvalue = rep(NA,ncol(cluster))
#
#
#   for(cli  in seq_len(max(cluster[1,]))){
#     sample <- which(cluster[1,]==cli)
#     stats <- distribution[1,sample]
#     pvalue[sample] <- compute_troendle(rbind(stats,distr_head[,seq_along(sample),drop=F]),
#                                        alternative = "greater")$main[,2]
#
#
#   }
#
#   list(main =cbind(statistic = distribution[1,], pvalue = pvalue))
#
#
# }
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


