#' @export
print.lmpermutation_table<-function(x, digits = 4, na.print = "", ...){
  cat(attr(x,"heading"), sep = "\n\n")
  cat(attr(x,"type"), sep = "\n\n")
  if(is.data.frame(x)){
    print(as.matrix(x), digits = digits, na.print = na.print, ...)
  }else{
    print.listof(x, digits = digits, ...)
  }
}