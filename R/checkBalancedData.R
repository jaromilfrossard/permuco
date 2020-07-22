#' @importFrom stats formula
checkBalancedData <- function(fixed_formula,data){
  tfixed=terms(fixed_formula,data=data)
  ho=attr(tfixed,"term.labels")[attr(tfixed,"order")==max(attr(tfixed,"order"))][1]
  mmho=model.matrix(formula(paste("~",ho,collapse="")),data=data)
  dataPerGroup=colSums(matrix(mmho[,attr(mmho,"assign")==1],nrow=NROW(mmho)))
  if(!(sum(!(dataPerGroup==dataPerGroup[1]))==0)){
    warning("The data are not balanced, the results may not be exact.")
  }
}