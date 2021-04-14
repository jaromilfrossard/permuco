############################################################################################################
#' Cluster-depth correction
#'
#' @description Compute the clusterdepth test correction given a matrix a permuted statistical signals.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param threshold A scalar that represents the threshold to create the clusters.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @param depth_scale A character string indicating the measure of the depth. Default is \code{"head_and_tail"} corrects independently from head and tail. \code{"unique_depth"} is also available and use a unique depth measure.
#' @param border A character string indicating the method to handle the clusters at the border. Default is \code{"reverse"} and will reverse the measure of the depth. \code{"ignore"} is also available and consider these cluster similarly to the others.
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





  pvalue <- rep(NA,ncol(cluster))


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




