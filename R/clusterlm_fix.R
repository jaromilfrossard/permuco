#' @importFrom stats rnorm qf qt pt pf
clusterlm_fix <- function(formula, data, method, test, threshold, np, P, rnd_rotation, aggr_FUN, E, H,
                          cl, multcomp, alpha, p_scale, coding_sum, ndh, return_distribution, new_method){


  multcomp <- match.arg(multcomp, c("clustermass", "tfce", "troendle","minP", "bonferroni", "holm", "benjamini_hochberg"),
                        several.ok = T)


  if(is.null(method)){method <- "freedman_lane"}

  if(!new_method){
    method <- match.arg(method,c("freedman_lane","kennedy","huh_jhun","manly",
                                "terBraak","draper_stoneman","dekker"))
  }

  if(is.null(aggr_FUN)){
    switch(test,
           "t"={
             fun_name <- "the sum of squares"
             aggr_FUN <- function(x)sum(x^2)},
           "fisher"={
             fun_name <- "the sum"
             aggr_FUN <- function(x)sum(x)})
  }else{
    fun_name <- "a user-defined function"
  }




  switch(method,
    "draper_stoneman" = {
      funP <- function(...) {
        cluster_draper_stoneman(...)
      }},
    "manly" = {
      funP <- function(...) {
        cluster_manly(...)
      }},
    "huh_jhun" = {
      funP <- function(...) {
        cluster_huh_jhun(...)
      }},
    "freedman_lane" = {
      funP <- function(...) {
        cluster_freedman_lane(...)
      }},
    "kennedy" = {
      funP <- function(...) {
        cluster_kennedy(...)
      }},
    "dekker" = {
    funP <- function(...) {
      cluster_dekker(...)
    }},
    "terBraak" = {
    funP <- function(...) {
      cluster_terBraak(...)
    }},
    {warning(paste("the method",method, "is not defined. Choose between freedman_lane, huh_jhun, dekker, terBraak or see help."))
      funP <- function(...){eval(parse(text=paste("cluster_",test,"_",method,"(...)",sep="",collpase="")))}})

  #preprocess formula=========================================


  if(!(class(formula[[2]]) == "matrix")){
    formula[[2]] <- call("as.matrix", formula[[2]])}


  #reshape data
  #all
  #========
  m <- match(c("formula", "data"), names(cl), 0L)
  mf <- cl[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- as.formula(paste(c(formula[[2]],"~",formula[[3]]),collapse=""))
  mf[[1L]] <- quote(stats::model.frame)
  #========



  mf <- eval(mf, parent.frame(n=2))
  if(coding_sum){mf <- changeContrast(mf, contr = contr.sum)}


  mm <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  terms <- terms(formula, special = "Error", data = data)
  indError <- attr(terms, "specials")$Error

  ####select test==============================
  switch(test,
         "fisher" = {
           col_ref <- attr(mm,"assign")
           colx <- 1:max(attr(mm,"assign"))
           names(colx) <- attr(terms,"term.labels")
         },
         "t" = {
           col_ref <- 1:length(attr(mm,"assign"))
           qr_mm <- qr(mm)
           colx <- qr_mm$pivot[1:qr_mm$rank]
           if(method != "huh_jhun"){
             colx <- colx[(attr(mm,"assign")[colx])!=0]}
           names(colx) <- colnames(mm)[colx]})

  #permutation matrix====================================

  #check dim of P
  if (!is.null(P)) {
    check_P <- check_P(
      P = P,
      method = method,
      test = test,
      n = NROW(y),
      ncol_x2 = as.numeric(table(attr(mm,"assign")[attr(mm,"assign")!=0])),
      ncol_x = NCOL(mm)
    )
    if (!check_P) {
      np <- np(P)
      P <- NULL
      warnings("P argument is not valid and will be recomputed")
    }
  }

  #create permutation matrices ==============================
  if(is.null(P)){switch(method,
                        "huh_jhun" = {
                          switch (test,
                            "t" = {P <- Pmat(np = np, n = NROW(y) - NCOL(mm) + 1)},
                            "fisher" = {{
                              P <- lapply(as.numeric(table(col_ref))[-1],
                                          function(cx){Pmat(np = np, n = NROW(y) - NCOL(mm) + cx)})}}
                            )},
                        {P <- Pmat(np = np, n = NROW(y))})}


  if(sum(np(P) <= 1999)>0){
    warning("The number of permutations is below 2000, p-values might be unreliable.")
  }
  np <- np(P)

  #create rnd_rotation matrices==============================
  if(method=="huh_jhun" & is.null(rnd_rotation)){
    rnd_rotation <- matrix(rnorm(NROW(y)^2),ncol=NROW(y))
  }


  ###initialisze output
  multiple_comparison <- list()
  length(multiple_comparison) <- length(colx)
  names(multiple_comparison) <-names(colx)

    if(test == "t"){
    multiple_comparison_less <- multiple_comparison_greater <- list()
    length(multiple_comparison_less) <- length(multiple_comparison_greater) <- length(colx)
    names(multiple_comparison_less) <- names(multiple_comparison_greater) <-names(colx)
    }else{
      multiple_comparison_less <- multiple_comparison_greater <- NULL
    }

  ##adjust multiple threshold
  switch(test,
         "t" = {
           df <- compute_degree_freedom_fix(test = test,mm = mm,assigni = colx)},
         "fisher" = {
           df <- compute_degree_freedom_fix(test = test,mm = mm,assigni = attr(mm,"assign"))})


  if(is.null(threshold)){
    switch(test,
           "t" = {
             threshold <- qt(p = 0.975,df = df)},
           "fisher" = {
             threshold <- qf(p = 0.95, df1 = df[,1],df2 = df[,2])})
  }else if(length(threshold)==1){
    threshold <- rep(threshold,length(colx))
  }else if(length(threshold)>1){
    threshold <- as.numeric(matrix(threshold, nrow = length(colx)))
    }





  args <- list(y = y, mm = mm, P = P, rnd_rotation = rnd_rotation, test = test)



  for(i in 1:length(colx)){
    ##compute distribution
    args$colx <- which(col_ref == colx[i])
    switch(test,
      "t" = {dfi <- df[i]},
      "fisher" = {dfi <- df[i,]}
    )

    if(method == "huh_jhun"&test =="fisher"){args$P <- P[[i]]}
    distribution <- t(funP(args = args))
    pvalue <- apply(distribution,2,function(col)compute_pvalue(distribution = col))

    switch (test,
      "t" = {pvalue_para <- 2*pt(abs(distribution[1,]),df = dfi,lower.tail = F)},
      "fisher" = {pvalue_para <- pf(distribution[1,],df1 =  dfi[1],df2 =  dfi[2],lower.tail = F)})

    #####uncorrected
    multiple_comparison[[i]]$uncorrected <- list(main = cbind(statistic = distribution[1,],pvalue = pvalue,pvalue_para = pvalue_para),
                                                test_info = list(test = test, df = dfi, alternative = "two.sided", method = method, np = np,
                                                                 nDV = ncol(y), fun_name = fun_name))
    if(return_distribution){multiple_comparison[[i]]$uncorrected$distribution = distribution}

    ##pscale change
    if(p_scale){
      distribution0 <- distribution
      distribution <- distribution_to_pscale(distribution0, test = test, alternative = "two.sided")}

    #compute multiple comparison for two.sided
    multiple_comparison[[i]] <- c(multiple_comparison[[i]],
    switch_multcomp(multcomp = multcomp, distribution = distribution, threshold = threshold[i],aggr_FUN = aggr_FUN,
                    alternative = "two.sided", E = E,H = H,ndh =ndh,pvalue = pvalue, alpha = alpha))
    ##one-sided test
    if(test == "t"){
      #greater
      ##pscale change
      alternative <- "greater"
      if(p_scale){
        distribution <- distribution_to_pscale(distribution0, test = test, alternative = alternative)
        pvalue_para = pt(distribution0[1,],df = df[i],lower.tail = F)}else{
        pvalue_para = pt(distribution[1,],df = df[i],lower.tail = F)}

      pvalue <- apply(distribution,2,function(col)compute_pvalue(distribution = col,alternative = alternative))

      multiple_comparison_greater[[i]]$uncorrected <- list(main = cbind(statistic = distribution[1,],pvalue = pvalue,pvalue_para = pvalue_para),
                                                          test_info = list(test = test, df = df, alternative = alternative, method = method, np = np,
                                                                           nDV = ncol(y), fun_name = fun_name))
      multiple_comparison_greater[[i]] <- c(multiple_comparison_greater[[i]],
                                         switch_multcomp(multcomp = c("clustermass",multcomp[!multcomp%in%"tfce"]), distribution = distribution,
                                                         threshold = threshold[i],aggr_FUN = aggr_FUN,alternative = alternative,
                                                         E = E,H = H,ndh =ndh,pvalue = pvalue, alpha = alpha))
      #less
      ##pscale change
      alternative <- "less"
      if(p_scale){
        distribution <- distribution_to_pscale(distribution0, test = test, alternative = "less")
        alternative <- "greater"
        pvalue_para <- pt(distribution0[1,],df = dfi,lower.tail = T)}else{
          pvalue_para <- pt(distribution[1,],df = dfi,lower.tail = T)}

      pvalue <- apply(distribution,2,function(col)compute_pvalue(distribution = col,alternative = alternative))


      multiple_comparison_less[[i]]$uncorrected <- list(main = cbind(statistic = distribution[1,], pvalue = pvalue, pvalue_para = pvalue_para),
                                                       test_info = list(test = test, df = dfi, alternative = alternative, method = method, np = np,
                                                                        nDV = ncol(y), fun_name = fun_name))
      multiple_comparison_less[[i]] <- c(multiple_comparison_less[[i]],
                                         switch_multcomp(multcomp = c("clustermass",multcomp[!multcomp%in%"tfce"]),distribution = distribution,
                                                         threshold = threshold[i],aggr_FUN = aggr_FUN,alternative = alternative,
                                                         E = E, H = H, ndh =ndh, pvalue = pvalue, alpha = alpha))}
  }

  ####create table
  # cluster_table <- NULL
  # cluster_table_greater <- cluster_table_less <- NULL
  # cluster_table <- cluster_table(multiple_comparison)
  # cluster_table_greater <- cluster_table_less <- NULL
  # if(test=="t"){
  #   cluster_table_greater <- cluster_table(multiple_comparison_greater)
  #   cluster_table_less <- cluster_table(multiple_comparison_less)}

  ###parametic model
  mod_lm <- lm.fit(x = mm, y = y)

  #reshape args
  # multcomp = unique(c("uncorrected",
  #                     "clustermass",multcomp))[unique(c("uncorrected","clustermass",
  #                                                            multcomp))%in%c("uncorrected","clustermass",
  #                                                                            "tfce","troendle","bonferroni","holm",
  #                                                                            "benjamini_hochberg")]

  #output
  out <- list()
  out$y <- y
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
  out$test <- test
  out$threshold <- threshold
  out$P <- P
  out$np <- np
  out$df <- df
  out$rnd_rotation <- rnd_rotation
  out$multcomp <- multcomp
  out$multiple_comparison <- multiple_comparison
  out$multiple_comparison_greater <- multiple_comparison_greater
  out$multiple_comparison_less <- multiple_comparison_less
  # out$cluster_table = cluster_table
  # out$cluster_table_greater = cluster_table_greater
  # out$cluster_table_less = cluster_table_less
  out$alpha <- alpha
  out$method <- method
  out$fun_name <- fun_name
  class(out) <- "clusterlm"
  return(out)
}
