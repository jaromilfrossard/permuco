full_table = function(x, multcomp = NULL, ...){
  if(is.null(multcomp)){multcomp="clustermass"}
  dotargs = list(...)
  ct = lapply(1:length(x), function(j){
    effect = x[[j]]
    info = effect$uncorrected$test_info


    tab = effect[[multcomp]]$main[,c(1,2)]


    tab = data.frame(tab)

    colnames(tab) = c(info$test ,paste0("P(>) ",multcomp))

    attr(tab,"effect_name") = names(x)[j]
    attr(tab,"multcomp") = multcomp
    attr(tab,"method") = info$method
    attr(tab,"test") = info$test
    attr(tab,"alternative") = info$alternative
    attr(tab,"df") = info$df
    attr(tab,"np") = info$np
    class(tab) = append("cluster_table",class(tab))
    tab
  })
  class(ct) = append("listof_cluster_table",class(ct))
  names(ct) = names(x)
  ct



}