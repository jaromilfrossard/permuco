cluster_table = function(x, multcomp = NULL, alpha = 0.05, ...){

  if(is.null(multcomp)){multcomp="clustermass"}

  if(multcomp=="clustermass"){
    return(cluster_table_clustermass(x, ...))
  }



  dotargs = list(...)
  ct = lapply(1:length(x), function(j){
    effect = x[[j]]
    info = effect$uncorrected$test_info
    selected= effect[[multcomp]]$main[,2]< alpha
    cl = (selected-c(0,selected[-length(selected)]))!=0
    cl = cumsum(cl)+1
    unique_cluster = unique(cl)


    tab = t(sapply(unique_cluster,function(i){
      cl_select = cl == i
      timepoint = c(1:length(cl_select))[cl_select]
      c(range(timepoint),
        effect[[multcomp]]$main[timepoint[1],2]<alpha)
    }))

    tab = data.frame(tab)

    colnames(tab) = c("start","end", "P(>)")
    tab$`P(>)` = as.factor(tab$`P(>)`)
    levels(tab$`P(>)`)[levels(tab$`P(>)`)=="0"] = "n.s."
    levels(tab$`P(>)`)[levels(tab$`P(>)`)=="1"] = "sign"

    rownames(tab) = unique_cluster
    attr(tab,"effect_name") = names(x)[j]
    attr(tab,"multcomp") = multcomp
    attr(tab,"nDV") = info$nDV
    attr(tab,"method") = info$method
    attr(tab,"test") = info$test
    attr(tab,"alternative") = info$alternative
    attr(tab,"type") = info$type
    attr(tab,"df") = info$df
    attr(tab,"np") = info$np
    attr(tab,"table_type") = "cluster"
    class(tab) = append("multcomp_table",class(tab))
    tab
  })
  class(ct) = append("listof_multcomp_table",class(ct))
  names(ct) = names(x)
  ct
}