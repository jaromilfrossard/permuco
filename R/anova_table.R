# Anova Table
#
# @description Create an anova table.
# @param model A model of class lm
# @return A table with the sums of squares, the F statistics and p-values.
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom stats pf
anova_table <- function(model){
  switch(class(model),
         "lm"={
           mf <- model.frame(model)
           mm <- model.matrix(attr(mf, "terms"), data=mf)
           Y <- model.response(mf)
           assign=model$assign
           names <- attr(model$terms, "term.label")
           rss <- sum(model$residuals^2)
           df.residual=model$df.residual},
         "list"={
           mm<-qr.X(model$qr)
           assign<-model$assign
           res<-model$residuals
           Y<-model$fitted.values+res
           rss <- sum(res^2)
           names <- attr(model, "term.label")
           df.residual<-model$df.residual
           names <- attr(model$terms, "term.label")
         })


  table = t(sapply(1:max(assign),function(a){
    col = assign==a
    Df = sum(col)
    x1 = matrix(mm[,!col],ncol=sum(!col))
    x2=matrix(mm[,col],ncol=sum(col))
    qr_x1 <- qr(x1)
    qr_r1x2 <- qr(qr.resid(qr_x1,x2))
    Sum_Square=sum(qr.fitted(qr_r1x2,Y)^2)
    F_stat=(Sum_Square)/rss/Df*df.residual
    p_value=1-pf(F_stat,df1=Df,df2=df.residual)
    return(c(Sum_Square, Df, F_stat, p_value))}))
  table=rbind(table,c(rss,df.residual,NA,NA))
  table=as.data.frame(table)
  colnames(table)=c("SS","df","F","parametric P(>F)")
  row.names(table)=c(names,"Residuals")
  attr(table, "heading")= "Anova Table"
  attr(table, "type")= "parametric p-values"
  class(table)=c("lmpermutation_table","data.frame")
  return(table)
}