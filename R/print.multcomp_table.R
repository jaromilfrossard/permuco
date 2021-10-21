#' @export
print.multcomp_table <- function(x, ...) {
  cat("Effect: ",attr(x,"effect_name"), ".\n",sep="")
  cat("Alternative Hypothesis: ",attr(x,"alternative"), ".\n",sep="")
  cat("Statistic: ",attr(x,"test"),"(",paste(attr(x,"df"),collapse=", "),")", ".\n",sep="")
  cat("Resampling Method: ",attr(x,"method"), ".\n",sep="")
  cat("Type of Resampling: ",attr(x,"type"), ".\n",sep="")
  cat("Number of Dependant Variables: ",attr(x,"nDV"), ".\n",sep="")
  cat("Number of Resamples: ",attr(x,"np"), ".\n",sep="")
  cat("Multiple Comparisons Procedure: ",attr(x,"multcomp"), ".\n",sep="")
  if(attr(x,"multcomp") == "clustermass"){
    cat("Threshold: ",attr(x,"threshold"),".\n",sep="")
    cat("Mass Function: ",attr(x,"fun_name"),".\n",sep="")
  }
  if(attr(x,"table_type")=="cluster"){
    if(attr(x,"multcomp")=="clustermass"){
      cat("Table of clusters.\n")
    }else if(attr(x,"multcomp")!="clustermass"){
      cat("Table of pseudo-clusters.\n")
    }
  }else if(attr(x,"table_type")=="full"){
    cat("Table of all tests.\n")
  }

  cat("\n")
  if(!is.null(attr(x,"nocluster"))){
    if(attr(x,"nocluster")){
      cat("No cluster above the threshold.\n")}else{
        print.data.frame(x)}
  }else{
    print.data.frame(x)}
  cat("\n\n")
}