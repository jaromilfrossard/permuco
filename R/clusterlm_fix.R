#' @importFrom stats rnorm qf qt
clusterlm_fix <- function(formula, data, method, test, threshold, np, P, rnd_rotation, aggr_FUN, E, H,
                          cl, multcomp, alpha, p_scale, coding_sum, ndh, return_distribution, new_method){

  if(is.null(method)){method = "freedman_lane"}

  if(!new_method){
    method = match.arg(method,c("freedman_lane","kennedy","huh_jhun","manly",
                                "terBraak","draper_stoneman","dekker"))
  }



  switch(method,
    "draper_stoneman" = {
      funP = function(...) {
        cluster_draper_stoneman(...)
      }},
    "manly" = {
      funP = function(...) {
        cluster_manly(...)
      }},
    "huh_jhun" = {
      funP = function(...) {
        cluster_huh_jhun(...)
      }},
    "freedman_lane" = {
      funP = function(...) {
        cluster_freedman_lane(...)
      }},
    "kennedy" = {
      funP = function(...) {
        cluster_kennedy(...)
      }},
    "dekker" = {
    funP = function(...) {
      cluster_dekker(...)
    }},
    "terBraak" = {
    funP = function(...) {
      cluster_terBraak(...)
    }},
    {warning(paste("the method",method, "is not defined. Choose between freedman_lane, huh_jhun, dekker, terBraak or see help."))
      funP=function(...){eval(parse(text=paste("cluster_",test,"_",method,"(...)",sep="",collpase="")))}})

  #preprocess formula=========================================


  if(!(class(formula[[2]]) == "matrix")){
    formula[[2]] <- call("as.matrix", formula[[2]])}


  #reshape data
  #all
  #========
  m <- match(c("formula", "data"), names(cl), 0L)
  mf <- cl[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula = as.formula(paste(c(formula[[2]],"~",formula[[3]]),collapse=""))
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
           colx <- 1:length(attr(mm,"assign"))
           if(method != "huh_jhun"){
             colx <- colx[attr(mm,"assign")!=0]}
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
      np = np(P)
      P = NULL
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
                        {P = Pmat(np = np, n = NROW(y))})}


  if(sum(np(P) <= 1999)>0){
    warning("The number of permutations is below 2000, p-values might be unreliable")
  }
  np <- np(P)

  #create rnd_rotation matrices==============================
  if(method=="huh_jhun" & is.null(rnd_rotation)){
    rnd_rotation <- matrix(rnorm(NROW(y)^2),ncol=NROW(y))
  }


  ###initialisze output
  multiple_comparison <- list()
  length(multiple_comparison)<- length(colx)
  names(multiple_comparison) <-names(colx)

    if(test == "t"){
    multiple_comparison_left <- multiple_comparison_right <- list()
    length(multiple_comparison_left)<- length(multiple_comparison_right) <- length(colx)
    names(multiple_comparison_left)<- names(multiple_comparison_right) <-names(colx)
    }else{
      multiple_comparison_left <- multiple_comparison_right <- NULL
    }

  ##adjust multiple threshold
  if(is.null(threshold)){
    switch(test,
           "t" = {
             df = compute_degree_freedom_fix(test = test,mm = mm,assigni = colx)
             threshold = qt(p = 0.975,df = df)},
           "fisher" = {
             df = compute_degree_freedom_fix(test = test,mm = mm,assigni = attr(mm,"assign"))
             threshold = qf(p = 0.95, df1 = df[,1],df2 =df[,2])})
  }else if(length(threshold)==1){
    threshold = rep(threshold,length(colx))
  }else if(length(threshold)>1){
    threshold = as.numeric(matrix(threshold,nrow=length(colx)))
    }




  args <- list(y = y, mm = mm, P = P, rnd_rotation = rnd_rotation, test = test)

  for(i in 1:length(colx)){
    ##compute distribution
    args$colx <- which(col_ref == colx[i])
    if(method == "huh_jhun"&test =="fisher"){args$P = P[[i]]}
    distribution<- t(funP(args = args))
    pvalue <- apply(distribution,2,function(col)compute_pvalue(distribution = col))

    #####uncorrected
    multiple_comparison[[i]]$uncorrected = list(main = cbind(statistic = distribution[1,],pvalue = pvalue))
    if(return_distribution){multiple_comparison[[i]]$uncorrected$distribution = distribution}

    ##pscale change
    if(p_scale){
      distribution0 <- distribution
      distribution <- distribution_to_pscale(distribution0, test = test, lateraltiy = "bilateral")}

    #compute multiple comparison for bilateral
    multiple_comparison[[i]] = c(multiple_comparison[[i]],
    switch_multcomp(multcomp = c("maris_oostenveld",multcomp),distribution = distribution, threshold = threshold[i],aggr_FUN = aggr_FUN,
                    laterality = "bilateral", E = E,H = H,ndh =ndh,pvalue = pvalue, alpha = alpha))
    ##unilateral test
    if(test == "t"){
      #right
      ##pscale change
      lateraltiy = "right"
      if(p_scale){
        distribution <- distribution_to_pscale(distribution0, test = test, lateraltiy = lateraltiy)}

      pvalue <- apply(distribution,2,function(col)compute_pvalue(distribution = col,laterality = lateraltiy))
      multiple_comparison_right[[i]]$uncorrected = list(main = cbind(statistic = distribution[1,],pvalue = pvalue))
      multiple_comparison_right[[i]] = c(multiple_comparison_right[[i]],
                                         switch_multcomp(multcomp = c("maris_oostenveld",multcomp[!multcomp%in%"tfce"]), distribution = distribution,
                                                         threshold = threshold[i],aggr_FUN = aggr_FUN,laterality = lateraltiy,
                                                         E = E,H = H,ndh =ndh,pvalue = pvalue, alpha = alpha))
      #left
      ##pscale change
      lateraltiy = "left"
      if(p_scale){
        distribution <- distribution_to_pscale(distribution0, test = test, lateraltiy = "left")
        lateraltiy = "right"
      }
      pvalue <- apply(distribution,2,function(col)compute_pvalue(distribution = col,laterality = lateraltiy))
      multiple_comparison_left[[i]]$uncorrected = list(main = cbind(statistic = distribution[1,],pvalue = pvalue))
      multiple_comparison_left[[i]] = c(multiple_comparison_left[[i]],
                                         switch_multcomp(multcomp = c("maris_oostenveld",multcomp[!multcomp%in%"tfce"]),distribution = distribution,
                                                         threshold = threshold[i],aggr_FUN = aggr_FUN,laterality = lateraltiy,
                                                         E = E, H = H, ndh =ndh, pvalue = pvalue, alpha = alpha))}
  }

  ####create table
  cluster_table <- cluster_table(multiple_comparison)
  cluster_table_right <- cluster_table_left <- NULL
  if(test=="t"){
    cluster_table_right <- cluster_table(multiple_comparison_right)
    cluster_table_left <- cluster_table(multiple_comparison_left)}

  ###parametic model
  mod_lm = lm.fit(x =mm, y = y)

  #reshape args
  multcomp = unique(c("uncorrected",
                      "maris_oostenveld",multcomp))[unique(c("uncorrected","maris_oostenveld",
                                                             multcomp))%in%c("uncorrected","maris_oostenveld",
                                                                             "tfce","troendle","bonferroni","holm",
                                                                             "benjaminin_hochberg")]

  #output
  out = list()
  out$coefficients = mod_lm$coefficients
  out$residuals = mod_lm$residuals
  out$effects = mod_lm$effects
  out$rank = mod_lm$rank
  out$fitted.values = mod_lm$fitted.values
  out$assign = mod_lm$assign
  out$qr = mod_lm$qr
  out$df.residual = mod_lm$df.residual
  out$xlevels = mod_lm$xlevels
  out$terms = mod_lm$terms
  out$model = mod_lm$model
  out$model.matrix = mm
  out$test = test
  out$threshold = threshold
  out$P = P
  out$rnd_rotation = rnd_rotation
  out$multcomp = multcomp
  out$multiple_comparison = multiple_comparison
  out$multiple_comparison_right = multiple_comparison_right
  out$multiple_comparison_left = multiple_comparison_left
  out$cluster_table = cluster_table
  out$cluster_table_right = cluster_table_right
  out$cluster_table_left = cluster_table_left
  out$alpha = alpha
  out$method = method
  class(out) <- "clusterlm"
  return(out)
}
