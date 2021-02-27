#' @export
print.listof_multcomp_table<- function(x,...){
  dotargs = list(...)
  ei = NULL
  if(!is.null(dotargs$effect)){
    ei = which(names(x)%in%dotargs$effect)}
  if(length(ei)>0){
    for(i in ei){print(x[[i]])}
  }else{
    for(i in 1:length(x)){
      print(x[[i]])}}
}