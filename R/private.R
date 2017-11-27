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


#compute
anova_table_rnd = function(args){
  table<-t(sapply(1:max(attr(args$mm,"assign")),function(i){
    args$i = i
    effect_anova_rnd(args = args)}))
  table = as.data.frame(table)
  class(table) = c("lmpermutation_table","data.frame")
  return(table)
  }


effect_anova_rnd = function(args){
  #select x
  mm = args$mm
  assign = attr(mm,"assign")
  select_x = assign==args$i
  select_within = assign == (args$link[3,args$i])

  ##zmat
  z = khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
  z = qr.resid(qr(mm),z)
  qr_z = qr(z)

  qr_d = qr(mm[,!select_x,drop=F])
  rdx = qr.resid(qr_d,mm[,select_x,drop=F])


  qr_rdx = qr(rdx)
  qr_rdz = qr(qr.resid(qr_d,z))

  SSn = sum(qr.fitted(qr_rdx,args$y)^2)
  SSd = sum(qr.fitted(qr_rdz,args$y)^2)

  MSEn = SSn/ qr_rdx$rank
  MSEd = SSd/ qr_rdz$rank

  Fstat = MSEn/MSEd

  out= c(SSn = SSn, dfn = qr_rdx$rank,
           SSd = SSd, dfd = qr_rdz$rank,
           MSEn = MSEn, MSEd = MSEd,
         'F' = Fstat, 'parametric P(>F)' = 1-pf(q =Fstat, df1 = qr_rdx$rank, df2 =qr_rdz$rank))
  out
}



link=function(formula_f,formula_within){
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


#'@importFrom Matrix KhatriRao
khatrirao= function(a,b,c = NULL){
  if(!is.null(c)){
    if(dim(c)[2]==0){c=NULL}}
  if(dim(b)[2]==0){return(a)}

  out = t(as.matrix(KhatriRao(t(a),t(b))))
  if(!is.null(c)){
    out= t(as.matrix(KhatriRao(t(out),t(c))))
  }
  out}



#'@importFrom stats contrasts<-
changeContrast <- function(data,contr){
  isfact <- sapply(data, function(c){is.factor(c)||is.character(c)})
  for(i in c( 1:ncol(data))[isfact]){
    data[ ,i] <- as.factor(data[ ,i])
    data[ ,i] <- droplevels(data[ ,i])
    contrasts(data[ ,i]) <- contr
  }
  return(data)
}



check_P <- function(P, method, test, n, ncol_x2, ncol_x) {
  check_format <- function(p){
    out <- F
    switch(class(p),
           "matrix" = {out <- !is.logical(tryCatch(as.Pmat(P),error = function(e){T}))},
           "Pmat" = {out <- T})
    return(out)}

  nrow_n <- method != "huh_jhun"
  check <- F
  switch(paste(nrow_n),
         "TRUE" = {if(check_format(P)& NROW(P) == n){check <- T}},  ###NOT huh jhun control class et dimension
         "FALSE" = { ####huh jhhun
           switch(test,
                  "t" = {if(check_format(P)& NROW(P) == n-ncol_x+1){check <- T}},
                  "fisher" = {
                    switch(class(P),
                           "list" = {
                             if(all(c(sapply(P,function(p){check_format(p)}),
                                      sapply(P,function(p){NROW(p)})== n - ncol_x + ncol_x2))){check <- T}})})})
  return(check)
}

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



compute_pvalue <- function(distribution, stat = distribution[1], laterality="bilateral", na.rm = T, digits = 10){
  if(na.rm){distribution = as.numeric(na.omit(distribution))}
  distribution=round(distribution,digits = digits)
  stat=round(stat, digits = digits)
  stat=matrix(stat,nrow = 1)
  switch(laterality,
         "bilateral" = {apply(stat,2, function(val)mean(abs(distribution) >= abs(val) , na.rm = na.rm))},
         "left" = {apply(stat,2, function(val)mean(distribution <= val , na.rm = na.rm))},
         "right" = {apply(stat,2, function(val)mean(distribution >= val, na.rm = na.rm))})
}

#' @importFrom stats na.omit
compute_all_pvalue <- function(distribution, laterality="bilateral", na.rm = T, digits = 10){
  if(na.rm){distribution = na.omit(distribution)}
  distribution=round(distribution,digits = digits)
  np = length(distribution)
  switch(laterality,
         "bilateral" = {distribution = - abs(distribution)},
         "left" = {},
         "right" = {distribution = - distribution})
  ceiling(rank(distribution))/np
}




check_distribution <- function(distribution, digits = 10, n_unique = 30){
  unique_d = apply(distribution,2,function(d){
    length(unique(round(d,digits = digits)))
  })
  if(sum(unique_d<n_unique)>0){
    warning(paste("the distribution of ",paste(colnames(distribution)[unique_d<n_unique],sep =", ",collapse = ", "), " may be discrete.", sep=" ",collapse = " "))
  }
}







