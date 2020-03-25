#'@importFrom stats lm.fit rnorm
aovperm_fix <- function( formula, data, method, type, np, coding_sum, P, rnd_rotation, new_method = NULL){

  if(is.null(coding_sum)){coding_sum = T}

  if(is.null(new_method)){new_method = F}

  if(is.null(method)){method = "freedman_lane"}



  if(!new_method){
    method = match.arg(method,c("freedman_lane","kennedy","huh_jhun","manly",
                                "terBraak","draper_stoneman","dekker"))
  }



  #select method==================================================
  switch(method,
         "freedman_lane"={funP=function(...){fisher_freedman_lane(...)}},
         "kennedy"={funP=function(...){fisher_kennedy(...)}},
         "huh_jhun"={funP=function(...){fisher_huh_jhun(...)}},
         "manly"={funP=function(...){fisher_manly(...)}},
         "terBraak"={funP=function(...){fisher_terBraak(...)}},
         "draper_stoneman"={funP=function(...){fisher_draper_stoneman(...)}},
         "dekker"={funP=function(...){fisher_dekker(...)}},
         {warning(paste("The method",method, "is not defined. Choose between freedman_lane, huh_jhun, dekker, terBraak or see help. Set new_method = TRUE to allow user-defined functions."))
           funP=function(...){eval(parse(text=paste("fisher_",method,"(...)",sep="",collpase="")))}})


  #preprocess data=========================================
  mf <- model.frame(formula = formula, data = data)
  if(coding_sum){mf <- changeContrast(mf, contr = contr.sum)} #change the contrasts of factors
  mm <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  name <- colnames(mm)


  #parametric ANOVA=========================================
  terms <- terms(formula)
  mod_lm <- lm.fit(y = y, x = mm)
  mod_lm$terms <- terms
  anova_table = anova_table(mod_lm)
  #name <- attr(terms(formula), "term.labels")
  name <- attr(terms(formula), "term.labels")




  #check P========================================
    if(!is.null(P)){
    check_P <- check_P(P = P, method = method,
                       test = "fisher", n = length(y),
                       ncol_x2 = as.numeric(table(attr(mm,"assign")[attr(mm,"assign")!=0])), ncol_x = NCOL(mm))
    if(!check_P){
      np = np(P)
      type = attr(P,"type")
      P = NULL
      warning("P argument is not valid and will be recomputed.")
    }
  }


  #create permutation matrices==============================

  if(is.null(P)){
    switch(method,
           "huh_jhun" = {
             P <- lapply(unique(attr(mm,"assign"))[-1], function(assigni){
               Pmat(np = np, n = length(y) - NCOL(mm) + sum(attr(mm,"assign")==assigni),type  = type)})},
           {P = Pmat(np = np, n = length(y), type = type)})
  }

  if(sum(np(P) <= 1999)>0){
    warning("The number of permutations is below 2000, p-values might be unreliable.")
  }
  type = attr(P,"type")
  np <- np(P)

  #create rnd_rotation matrices==============================
  if(method=="huh_jhun" & is.null(rnd_rotation)){
    rnd_rotation <- matrix(rnorm(length(y)^2),ncol=length(y))
  }


  ####select test==============================
  colx <- 1:max(attr(mm,"assign"))
  names(colx) = attr(terms,"term.labels")

  ###compute distribution==============================
  args <- list(y = y, mm = mm, P = P, rnd_rotation = rnd_rotation)


  distribution <- sapply(colx,function(i){
    args$colx <- which(attr(mm,"assign")==i)
    if(method=="huh_jhun"){args$P = args$P[[i]]}
    out = funP(args = args)
    if(method=="huh_jhun"){length(out) = max(np(P))}
    out
  })

  distribution = matrix(distribution,nrow=np)
  colnames(distribution) = names(colx)



  ##compute p value
  permutation_pvalue = apply(distribution,2,function(d){compute_pvalue(distribution = d,alternative="two.sided", na.rm = T)})
  check_distribution(distribution = distribution, digits = 10, n_unique = 200)



  #output
  table = anova_table
  table$pValue_Permutation=c(permutation_pvalue, NA)
  colnames(table)[4:5]=c("parametric P(>F)","Re-sample P(>F)")

  if(is.matrix(P)){type = attr(P,"type")}else if(is.list(P)){type = attr(P[[1]],"type")}


  attr(table,"type") = paste0("Resample test using ",method," to handle nuisance variables and ",paste(np,sep=", ",collapse = ", "), " ",type ,"s.")

  out=list()
  out$coefficients <- mod_lm$coefficients
  out$residuals <- mod_lm$residuals
  out$effects <- mod_lm$effects
  out$rank <- mod_lm$rank
  out$fitted.values <- mod_lm$fitted.values
  out$assign <- mod_lm$assign
  out$qr <- mod_lm$qr
  out$df.residual <- mod_lm$df.residual
  out$xlevels <- mod_lm$xlevels
  out$terms <- mod_lm$terms
  out$model <- mod_lm$model
  out$model.matrix <- mm
  out$table <- table
  out$distribution <- distribution
  out$P <- P
  out$rnd_rotation <- rnd_rotation
  out$coding_sum <- coding_sum
  out$np = np
  out$method = method

  class(out)<-"lmperm"
  return(out)
}


