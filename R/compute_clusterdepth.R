############################################################################################################
#' Cluster-depth correction
#'
#' @description Compute the clusterdepth test correction given a matrix a permuted statistical signals.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param threshold A scalar that represents the threshold to create the clusters.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @export
#' @family multcomp
compute_clusterdepth <- function(distribution,threshold, alternative = "two.sided"){
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))


  cluster <- get_cluster(distribution = distribution[1,,drop=F],
                         threshold = threshold, alternative = alternative,side= "all")

  cluster_head <- get_cluster(distribution = distribution, threshold = threshold, alternative = alternative, side = "starting")
  depth_head <- get_clusterdepth_head(cluster_head, border = "ignore")
  distr_head <- depth_distribution(distribution, head_mat = depth_head)


  pvalue_head <- rep(NA,ncol(cluster_head))
  #max_cl_size <- max(table(cluster_head[1,cluster_head[1,]!=0]))


  for(cli  in seq_len(max(cluster_head[1,]))){
    sample <- which(cluster_head[1,]==cli)
    di <- distr_head
    di <- rbind(c(distribution[1, sample],rep(0,ncol(di)-length(sample))),di)
    pvalue_head[sample] <- compute_troendle(di,
                                       alternative = alternative)$main[seq_along(sample),2]


  }


  cluster_tail <- get_cluster(distribution = distribution, threshold = threshold, alternative = alternative, side = "ending")
  depth_tail <- get_clusterdepth_tail(cluster_tail,border = "ignore")
  distr_tail <- depth_distribution(distribution, tail_mat = depth_tail)



  pvalue_tail <- rep(NA,ncol(cluster_tail))
  #max_cl_size <- max(table(cluster_tail[1,cluster_tail[1,]!=0]))


  for(cli  in seq_len(max(cluster_tail[1,]))){
    sample <- which(cluster_tail[1,]==cli)
    di <- distr_tail
    di <- rbind(c(rep(0,ncol(di)-length(sample)),distribution[1, sample]),di)
    pvalue_tail[sample] <- compute_troendle(di,
                                       alternative = alternative)$main[seq_along(sample)+ncol(distr_tail)-length(sample),2]

    # stats <- distribution[1,sample]
    # stats <- c(rep(0,max_cl_size-length(stats)),stats)
    # pvalue_tail[sample] <- compute_troendle(rbind(stats,distr_tail[,seq_len(max_cl_size),drop=F]),
    #                                         alternative = alternative)$main[seq_along(sample),2]


  }

  pvalue = apply(cbind(pvalue_head,pvalue_tail),1,max,na.rm = FALSE)

  list(main =cbind(statistic = distribution[1,], pvalue = pvalue, cluster_id = as.integer(cluster)),threshold = threshold)


}




