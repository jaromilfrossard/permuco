full_table = function(x, multcomp = NULL, ...){
  if(is.null(multcomp)){multcomp="clustermass"}
  dotargs = list(...)
  ct = lapply(1:length(x), function(j){
    effect = x[[j]]
    info = effect$uncorrected$test_info



    if(multcomp%in%c("clustermass","tfce")){
      tab = cbind(effect$uncorrected$main[,c(1)],
                  effect[[multcomp]]$main[,c(1,2)])
      tab = data.frame(tab)
      colnames(tab) = c(info$test,paste0("modified ",info$test) ,paste0("P(>) ",multcomp))
    }else{
      tab = cbind(effect$uncorrected$main[,c(1)],
                  effect[[multcomp]]$main[,c(2)])
      tab = data.frame(tab)
      colnames(tab) = c(info$test ,paste0("P(>) ",multcomp))
    }

    attr(tab,"effect_name") = names(x)[j]
    attr(tab,"multcomp") = multcomp
    attr(tab,"nDV") = info$nDV
    attr(tab,"method") = info$method
    attr(tab,"test") = info$test
    attr(tab,"alternative") = info$alternative
    attr(tab,"df") = info$df
    attr(tab,"np") = info$np
    attr(tab,"table_type") = "full"
    if(multcomp=="clustermass"){
      attr(tab,"threshold") = effect$clustermass$threshold
      attr(tab,"fun_name") = info$fun_name
    }
    class(tab) = append("multcomp_table",class(tab))
    tab
  })
  class(ct) = append("listof_multcomp_table",class(ct))
  names(ct) = names(x)
  ct



}