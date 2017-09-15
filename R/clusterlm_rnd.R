clusterlm_rnd <- function(formula, data, method, test, coding_sum, threshold, np, P, rnd_rotation,
                          aggr_FUN, E, H, cl, multcomp, alpha, p_scale, return_distribution, ndh, new_method){

  if(is.null(method)){method = "Rd_kheradPajouh_renaud"}

  if(!new_method){
    method = match.arg(method,c("Rd_kheradPajouh_renaud","Rde_kheradPajouh_renaud"))
  }


  switch(
    paste(method,sep="_"),
    "Rd_kheradPajouh_renaud" = {
      funP = function(...) {
        cluster_Rd_kheradPajouh_renaud_rnd(...)
      }},
    "Rde_kheradPajouh_renaud" = {
      funP = function(...) {
        cluster_Rde_kheradPajouh_renaud_rnd(...)
      }},
    {warning(paste("the method",method, "is not defined. Choose between Rd_kheradPajouh_renaud or Rde_kheradPajouh_renaud."))
      funP=function(...){eval(parse(text=paste("cluster_fisher_",method,"_rnd(...)",sep="",collpase="")))}})

  if(!(class(formula[[2]]) == "matrix")){
    formula[[2]] <- call("as.matrix", formula[[2]])}


  #Formula transforamtion
  terms<-terms(formula,special="Error",data=data)
  ind_error <- attr(terms, "specials")$Error
  error_term <- attr(terms, "variables")[[1 + ind_error]]
  formula_f <- update(formula, paste(". ~ .-",deparse(error_term, width.cutoff = 500L, backtick = TRUE)))
  e_term <- deparse(error_term[[2L]], width.cutoff = 500L,backtick = TRUE)

  #random/fix formula

  formula_allfixed <- as.formula(paste(c(formula_f[[2]],"~",formula_f[[3]],"+",e_term),collapse=""))
  formula_allfixed_design <- as.formula(paste(c("~",formula_f[[3]],"+",e_term),collapse=""))

  formula_within <- formula(paste("~", e_term, collapse = ""))
  formula_within<- formula(paste("~",deparse(error_term[[2]][[3]]),collapse=""))
  formula_id <- formula(paste("~",deparse(error_term[[2]][[2]]),collapse = ""))

  #Model frame
  mf <- model.frame(formula = formula_allfixed, data = data)
  mf_design <- model.frame(formula = formula_allfixed_design, data = data)
  if(coding_sum){mf_design <- changeContrast(mf_design, contr = contr.sum)}

  mf_f <- model.frame(formula = formula_f, data = mf_design)
  mf_id <- model.frame(formula = formula_id, data = as.data.frame(lapply(mf_design,function(col){
    col = as.factor(col)
    contrasts(col) = contr.sum
    col})))

  ##response
  mf <- eval(mf, parent.frame(n=2))
  y <- model.response(mf)

  ##link fixed random
  link = link(formula_f=formula_f,formula_within=formula_within)

  ###model .matrix
  mm_f <- model.matrix(attr(mf_f, "terms"), data = mf_f)
  mm_id <- model.matrix(attr(mf_id, "terms"), data = mf_id)[,-1,drop=F]
  name <- colnames(mm_f)

  ##checkk data


  checkBalancedData(fixed_formula = formula_f, data = cbind(mf))

  #compute permutation
  if (is.null(P)) {P = Pmat(np = np, n = NROW(y))}

  ##distribution
  args <- list(y = y, mm = mm_f, mm_id = mm_id, link = link, P = P)
  ###initialisze output
  multiple_comparison <- list()
  length(multiple_comparison) <- max(attr(mm_f,"assign"))
  names(multiple_comparison) <- attr(attr(mf_f,"terms"),"term.labels")

  for(i in 1:max(attr(mm_f,"assign"))){
    args$i = i
    distribution = funP(args = args)
    pvalue <- apply(distribution,2,function(col)compute_pvalue(distribution = col))

    #####uncorrected
    multiple_comparison[[i]]$uncorrected = list(main = cbind(statistic = distribution[1,],pvalue = pvalue))
    if(return_distribution){multiple_comparison[[i]]$uncorrected$distribution = distribution}

    ##pscale change
    if(p_scale){
      distribution0 <- distribution
      distribution <- distribution_to_pscale(distribution0, test = test, lateraltiy = "bilateral")
    }

    multiple_comparison[[i]] = c(multiple_comparison[[i]],
                                 switch_multcomp(multcomp = c("maris_oostenveld",multcomp),
                                                 distribution = distribution, threshold = threshold,
                                                 aggr_FUN = aggr_FUN, laterality = "bilateral", E = E,
                                                 H = H,ndh =ndh,pvalue = pvalue, alpha = alpha))}

  cluster_table <- cluster_table(multiple_comparison)

  ##sort effect
  cluster_table = cluster_table[order(link[3,], link[1,])]
  multiple_comparison = multiple_comparison[order(link[3,], link[1,])]

  attr(cluster_table,"type") <- paste("Permutation test using ",method," to handle noise variable and ",np, " permutations.")

  out=list()
  out$y = y
  out$model.matrix = mm_f
  out$model.matrix_id = mm_id = mm_id
  out$link = link
  out$P = P
  out$cluster_table = cluster_table
  out$multiple_comparison = multiple_comparison
  out$data=mf
  out$method = method
  out$alpha = alpha
  out$multcomp = multcomp
  out$threshold = threshold
  out$test = test
  class(out) <- "clusterlm"
  return(out)

}
