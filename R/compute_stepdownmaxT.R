#'  Step-down version of the max-T correction
#'
#' @description Compute the Step-down version of the max-T given a matrix a permuted statistics.
#'
#' @param distribution A matrix of permuted statistical signal. The first row indicating the observed statistics.
#' @param alternative A character string indicating the alternative hypothesis. Default is \code{"greater"}. Choose between \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#'
#' @export
#' @family multcomp
compute_stepdownmaxT = function (distribution, alternative) {
  alternative <- match.arg(alternative, c("two.sided","greater","less"))
  switch (alternative,
          "greater" = {
            test_rank <- rank(-distribution[1,])
            test_order <- sort(unique(test_rank))

            p_corrected <-
              lapply(seq_along(test_order),function(testi){
                ri <- test_order[testi]
                coltesti  <- which(test_rank==ri)
                coldistri <-which(test_rank>=ri)
                pvali <- sapply(distribution[1,coltesti],
                                function(stati){
                                  di = apply(distribution[,which(test_rank>=ri),drop=F],1,max)
                                  compute_pvalue(di, stat = stati,alternative= alternative)})
                cbind(cols = coltesti, pvalues = pvali)
              })

            },
          "less" = {
            test_rank <- rank(distribution[1,])
            test_order <- sort(unique(test_rank))

            p_corrected <-
              lapply(seq_along(test_order),function(testi){
                ri <- test_order[testi]
                coltesti  <- which(test_rank==ri)
                coldistri <-which(test_rank>=ri)
                pvali <- sapply(distribution[1,coltesti],
                                function(stati){
                                  di = apply(distribution[,which(test_rank>=ri),drop=F],1,min)
                                  compute_pvalue(di, stat = stati,alternative= alternative)})
                cbind(cols = coltesti, pvalues = pvali)
              })

          },
          "two.sided" = {
            test_rank <- rank(-abs(distribution[1,]))
            test_order <- sort(unique(test_rank))

            p_corrected <-
              lapply(seq_along(test_order),function(testi){
                ri <- test_order[testi]
                coltesti  <- which(test_rank==ri)
                coldistri <-which(test_rank>=ri)
                pvali <- sapply(distribution[1,coltesti],
                                function(stati){
                                  di = apply(abs(distribution[,which(test_rank>=ri),drop=F]),1,max)
                                  compute_pvalue(di, stat = stati,alternative= alternative)})
                cbind(cols = coltesti, pvalues = pvali)
              })


          }
  )

  p_corrected = do.call("rbind",p_corrected)

  cummax <- cummax(p_corrected[order(p_corrected[,1]),2])
  cummax = cummax[p_corrected[,1]]
  p_corrected = cbind(p_corrected,cummax = cummax)



  out <- list(main = cbind(statistic = distribution[1, ], pvalue = p_corrected[,3]))
  return(out)
}