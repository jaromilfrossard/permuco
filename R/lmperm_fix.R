#' @importFrom stats model.frame model.matrix model.response lm rnorm
lmperm_fix <- function(formula, data, method, np, P, rnd_rotation, new_method = NULL) {


  if(is.null(new_method)){new_method = F}
  if(is.null(method)){method = "freedman_lane"}
  if(!new_method){
    method = match.arg(method,c("freedman_lane","kennedy","huh_jhun","manly",
                                "terBraak","draper_stoneman","dekker"))
  }

  #select method==================================================
  switch(method,
         "freedman_lane"={funP=function(...){t_freedman_lane(...)}},
         "kennedy"={funP=function(...){t_kennedy(...)}},
         "huh_jhun"={funP=function(...){t_huh_jhun(...)}},
         "manly"={funP=function(...){t_manly(...)}},
         "terBraak"={funP=function(...){t_terBraak(...)}},
         "draper_stoneman"={funP=function(...){t_draper_stoneman(...)}},
         "dekker"={funP=function(...){t_dekker(...)}},
         {warning(paste("the method",method, "is not defined. Choose between freedman_lane, huh_jhun, dekker, terBraak or see help."))
           funP=function(...){eval(parse(text=paste("t_",method,"(...)",sep="",collpase="")))}})

  #preprocess data=========================================
  mf <- model.frame(formula=formula, data=data)
  mm <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  name <- colnames(mm)

  # #parametric LM=========================================
  mod_lm <- lm(formula = formula,data = data)
  table <- summary(mod_lm)$coefficients

  #check P argument
  if(!is.null(P)){
    check_P <- check_P(P = P, method = method,
                       test = "t", n = length(y),
                       ncol_x2 = rep(NCOL(mm)-1,NCOL(mm)), ncol_x = NCOL(mm))
    if(!check_P){
      np = np(P)
      P = NULL
      warnings("P argument is not valid and will be recomputed")
    }
  }

  #create permutation matrices==============================
  if(is.null(P)){switch(method,
                        "huh_jhun" = {
                          P <- Pmat(np = np, n = length(y) - NCOL(mm) + 1)},
                        {P = Pmat(np = np, n = length(y))})}
  if(sum(np(P) <= 1999)>0){
    warning("The number of permutations is below 2000, p-values might be unreliable")
  }

  ####create rnd_rotation matrices==============================
  if(method=="huh_jhun" & is.null(rnd_rotation)){
    rnd_rotation <- matrix(rnorm(length(y)^2),ncol=length(y))
  }

  ####select test==============================
  # qr decomposition take care of the rank < dim
  qr_mm = qr(mm)
  colx <- which(qr_mm$pivot<=qr_mm$rank)
  if(method != "huh_jhun"){colx <- colx[attr(mm,"assign")!=0]}
  names(colx) = colnames(mm)[colx]

  ###compute distribution==============================
  args <- list(y = y, mm = mm, P = P, rnd_rotation = rnd_rotation)

  distribution <- sapply(colx,function(i){
    args$colx <- i
    funP(args = args)
  })

  ### compute p_value============================
  r_pvalue <- l_pvalue <- bi_pvalue <- rep(NA, length(attr(mm,"assign")))

  r_pvalue[colx] = apply(distribution,2,function(d){compute_pvalue(distribution = d,laterality="right", na.rm = T)})
  l_pvalue[colx] = apply(distribution,2,function(d){compute_pvalue(distribution = d,laterality="left", na.rm = T)})
  bi_pvalue[colx] = apply(distribution,2,function(d){compute_pvalue(distribution = d,laterality="bilateral", na.rm = T)})

  check_distribution(distribution = distribution, digits = 10, n_unique = 200)

######create table===========================
  table <-
    as.data.frame(cbind(table, l_pvalue ,r_pvalue,bi_pvalue))
  colnames(table)[3:7] = c("t value","parametric Pr(>|t|)",
                           "permutation Pr(<t)","permutation Pr(>t)","permutation Pr(>|t|)")
  class(table) <- c("lmpermutation_table","data.frame")
  attr(table,"heading") <-
    c("Table of marginal t-test of the betas")
  attr(table,"type") = paste("Permutation test using",method,"to handle noise variable and",np, "permutations.")

#####output==================================
  out <- list()
  out$coefficients = mod_lm$coefficients
  out$residuals = mod_lm$residuals
  out$effects = mod_lm$effects
  out$rank = mod_lm$rank
  out$fitted.values = mod_lm$fitted.values
  out$assign = mod_lm$assign
  out$qr = mod_lm$qr
  out$df.residual = mod_lm$df.residual
  out$contrasts = mod_lm$contrasts
  out$xlevels = mod_lm$xlevels
  out$terms = mod_lm$terms
  out$model = mod_lm$model
  out$model.matrix = mm
  out$table = table
  out$distribution = distribution
  out$P = P
  out$rnd_rotation = rnd_rotation

  class(out) <- "lmperm"
  return(out)

}
