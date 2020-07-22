link <- function(formula_f,formula_within){
  t_f=terms(formula_f)
  t_within=terms(formula_within)
  factors_f = attr(t_f,"factors")[-1,,drop=F]
  factors_between=colSums(factors_f[rownames(factors_f)%in%attr(t_within,"term.labels"),,drop=F])==0
  factors_link=factors_f
  factors_link[!rownames(factors_f)%in%attr(t_within,"term.labels"),]=0
  factors_link = apply(factors_link,2,function(l){
    out=which(apply(factors_f,2,function(f)identical(as.logical(f),as.logical(l))))
    if(length(out)==0){out=0}
    out})
  link = rbind(order = attr(t_f,"order"),between_effect = factors_between*c(1:length(factors_between)),within = factors_link )
  link}