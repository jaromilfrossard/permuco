#' @importFrom stats aov update as.formula
aovperm_rnd <- function( formula, data, method, np, P, coding_sum, rnd_rotation, new_method = NULL){

  if(is.null(coding_sum)){coding_sum = T}

  if(is.null(new_method)){new_method = F}
  if(is.null(method)){method = "Rd_kheradPajouh_renaud"}
  if(!new_method){
    method = match.arg(method,c("Rd_kheradPajouh_renaud","Rde_kheradPajouh_renaud"))
  }



  switch(method,
         "Rd_kheradPajouh_renaud"={funP=function(...){fisher_Rd_kheradPajouh_renaud_rnd(...)}},
         "Rde_kheradPajouh_renaud"={funP=function(...){fisher_Rde_kheradPajouh_renaud_rnd(...)}},
         {warning(paste("the method",method, "is not defined. Choose between Rd_kheradPajouh_renaud or Rde_kheradPajouh_renaud. Set new_method = TRUE to allow user-defined functions."))
           funP=function(...){eval(parse(text=paste("fisher_",method,"_rnd(...)",sep="",collpase="")))}})

  #Formula transforamtion
  terms<-terms(formula,special="Error",data=data)
  ind_error <- attr(terms, "specials")$Error
  error_term <- attr(terms, "variables")[[1 + ind_error]]
  formula_f <- update(formula, paste(". ~ .-",deparse(error_term, width.cutoff = 500L, backtick = TRUE)))
  e_term <- deparse(error_term[[2L]], width.cutoff = 500L,backtick = TRUE)

  #random/fix formula

  formula_allfixed <- as.formula(paste(c(formula_f[[2]],"~",formula_f[[3]],"+",e_term),collapse=""))

  formula_within <- formula(paste("~", e_term, collapse = ""))
  formula_within<- formula(paste("~",deparse(error_term[[2]][[3]]),collapse=""))
  formula_id <- formula(paste("~",deparse(error_term[[2]][[2]]),collapse = ""))

  #Model frame
  mf <- model.frame(formula = formula_allfixed, data = data)
  if(coding_sum){mf <- changeContrast(mf, contr = contr.sum)}

  mf_f <- model.frame(formula = formula_f, data = mf)
  mf_id <- model.frame(formula = formula_id, data = as.data.frame(lapply(mf,function(col){
    col = as.factor(col)
    contrasts(col) = contr.sum
    col})))

  ##response
  y <- model.response(mf)

  ##link fixed random
  link = link(formula_f=formula_f,formula_within=formula_within)


  ###model .matrix
  mm_f <- model.matrix(attr(mf_f, "terms"), data = mf_f)
  mm_id <- model.matrix(attr(mf_id, "terms"), data = mf_id)[,-1,drop=F]

  name <- colnames(mm_f)

  ##checkk data
  checkBalancedData(fixed_formula = formula_f, data = cbind(y,mf))

  #compute permutation
  if (is.null(P)) {P = Pmat(np = np, n = length(y))}
  np = np(P)

  ##distribution
  args <- list(y = y, mm = mm_f, mm_id = mm_id, link = link, P = P)
  distribution<-sapply(1:max(attr(mm_f,"assign")),function(i){
    args$i = i
    funP(args = args)})
  distribution = matrix(distribution,nrow=np)
  colnames(distribution) = attr(attr(mf_f,"terms"),"term.labels")
  check_distribution(distribution = distribution, digits = 10, n_unique = 200)

  #anova table
  table = anova_table_rnd(args)
  rownames(table) = attr(attr(mf_f,"terms"),"term.labels")

  #permutation pvalue
  permutation_pvalue = apply(distribution,2,function(d){compute_pvalue(distribution = d,laterality="bilateral", na.rm = T)})
  table$'permutation P(>F)' = permutation_pvalue

  ##sort effect
  table = table[order(link[3,], link[1,]),]
  distribution = distribution[,order(link[3,], link[1,])]

  attr(table,"type") <- paste("Permutation test using",method,"to handle nuisance variables and",np, "permutations.")

  out=list()
  out$y = y
  out$model.matrix = mm_f
  out$model.matrix_id = mm_id = mm_id
  out$link = link
  out$P = P
  out$np = np
  out$table = table
  out$distribution = distribution
  out$data=mf
  out$method = method
  out$coding_sum = coding_sum
  class(out)<-"lmperm"
  return(out)

}
