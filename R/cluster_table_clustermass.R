
cluster_table_clustermass = function(x, ...){
  dotargs <- list(...)
  ct <- lapply(1:length(x), function(j){
    effect <- x[[j]]
    info <- effect$uncorrected$test_info
    unique_cluster <- unique(effect$clustermass$main[,3])
    unique_cluster <- unique_cluster[unique_cluster!=0]
    if(length(unique_cluster)==0){
      tab = data.frame()
      attr(tab,"nocluster") = T
    }else{
      tab <- t(sapply(unique_cluster,function(i){
        cl_select <- effect$clustermass$main[,3] == i
        timepoint <- c(1:length(cl_select))[cl_select]
        c(timepoint[1],timepoint[length(timepoint)],
          effect$clustermass$main[timepoint[1],1],
          effect$clustermass$main[timepoint[1],2])

      }))
      tab <- data.frame(tab)
      colnames(tab) <- c("start","end", "cluster mass", "P(>mass)")
      rownames(tab) <- unique_cluster
      attr(tab,"nocluster") <- F
    }
    attr(tab,"threshold") <- effect$clustermass$threshold
    attr(tab,"fun_name") <- info$fun_name
    attr(tab,"effect_name") <- names(x)[j]
    attr(tab,"multcomp") <- "clustermass"
    attr(tab,"nDV") <- info$nDV
    attr(tab,"method") <- info$method
    attr(tab,"test") <- info$test
    attr(tab,"type") = info$type
    attr(tab,"alternative") <- info$alternative
    attr(tab,"df") <- info$df
    attr(tab,"np") <- info$np
    attr(tab,"table_type") <- "cluster"

    class(tab) <- append("multcomp_table", class(tab))
    tab
  })
  class(ct) <- append("listof_multcomp_table", class(ct))
  names(ct) <- names(x)
  ct
}