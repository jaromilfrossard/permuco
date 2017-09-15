############################################################################################################
compute_maris_oostenveld <-  function(distribution, threshold, aggr_FUN, laterality = "right"){
    switch(laterality,
           "right" = {
             threshold <- abs(threshold)
             selected <- distribution > threshold
             extreme = function(x)max(x,na.rm = T)
           },
           "left" = {
             threshold <- -abs(threshold)
             selected <- distribution < threshold
             extreme = function(x)max(x,na.rm = T)},
           "bilateral" = {
             distribution <- abs(distribution)
             threshold <- abs(threshold)
             selected <- distribution > threshold
             extreme = function(x)max(x,na.rm = T)
           })

  ##create connected component labeling
  cl = (selected-cbind(0,selected[,-NCOL(selected)]))==1
  cl = t(apply(cl,1,cumsum))
  cl[!selected] = 0

  mass_distribution = sapply(1:NROW(cl),function(rowi){
    extreme(sapply(1:max(1,max(cl[rowi,])),function(i){
      aggr_FUN(distribution[rowi,cl[rowi,] == i])
    }))
  })
  mass_statistic =sapply(1:max(cl[1,]),function(i){
      aggr_FUN(distribution[1,cl[1,] == i])
    })
  statistic = c(NA,mass_statistic)[cl[1,]+1]


  pvalue = sapply(mass_statistic,function(mi)
    compute_pvalue(stat = mi,distribution = mass_distribution, laterality = "right"))
  pvalue = c(NA,pvalue)[cl[1,]+1]
  out = list(main = cbind(statistic = statistic, pvalue = pvalue, cluster_id = cl[1,]),
              distribution = mass_distribution)
  return(out)
}



############################################################################################################


compute_tfce_statistic <- function(signal, E, H, dh, dhi){
  ##reduce dhi to signal max
  dhi = dhi[dhi<=max(signal)]
  ##compute cluster by dh
  cl_log <- t(sapply(dhi,function(hi){signal > hi}))
  cl_beg <- (cl_log-cbind(F,cl_log[,-NCOL(cl_log), drop=F]))==1
  cl_num <- t(apply(cl_beg,1,cumsum));cl_num[!cl_log]<-0 ##unique integer by cluster by row
  extend_by_h=t(apply(cl_num,1,function(cl_row){
    c(0,sapply(1:max(cl_row),function(i){
      (sum(cl_row==i))^E
    }))[cl_row+1]
  }))
  ###multiple extend by height
  statistic = colSums(extend_by_h*dhi^H)*dh
  statistic
}


compute_tfce <- function(distribution, E, H, ndh){
  range_d <- range(abs(distribution))
  dhi <- seq(from = range_d[1], to =range_d[2], length.out = ndh)
  dh <- dhi[2]-dhi[1]

  distribution_tfce <- apply(abs(distribution),1, function(rowi){
    max(compute_tfce_statistic(signal = rowi, E = E, H = H,dh = dh, dhi = dhi))
  }
  )

  statistic <- compute_tfce_statistic(signal = abs(distribution[1,]), E = E, H = H,dh =dh, dhi = dhi)

  pvalue <- sapply(statistic,function(s)compute_pvalue(distribution = distribution_tfce,stat=s))

  out = list(main = cbind(statistic = statistic, pvalue = pvalue),
              distribution = distribution_tfce,
              dhi = dhi)
  return(out)



}


############################################################################################################

compute_bonferroni <- function(statistic = NULL,pvalue){
  pvalue <- pmin(pvalue*length(pvalue),1)
  out = list(main = cbind(statistic = statistic,pvalue = pvalue))
  return(out)
}

############################################################################################################

compute_holm <- function(statistic = NULL,pvalue){
  n = length(pvalue)
  i <- seq_len(n)
  o <- order(pvalue)
  ro <- order(o)
  pvalue <- pmin(1, cummax((n - i + 1L) * pvalue[o]))[ro]
  out=list(main = cbind(statistic = statistic,pvalue = pvalue))
  return(out)

}

############################################################################################################
compute_troendle <- function(distribution, pvalue, alpha, laterality){
  distribution_rank = apply(distribution,2,function(col){compute_all_pvalue(col,laterality = laterality)})

  p_corrected <- rep(1,length(pvalue))
  rank_uncorr <- rank(pvalue)
  #loop to recompute minimal value
  for(urank in sort(unique(rank_uncorr))){
    which_test <- which(urank==rank_uncorr)
    pvali <- distribution_rank[,which(urank<=rank_uncorr)]
    distr_min <- apply(pvali,1,min)
    p_corri <- compute_pvalue(distribution = distr_min, stat = matrix(pvalue[which_test],nrow=1), laterality = "left")
    p_corrected[which_test] <- p_corri
    if(sum(p_corri > alpha)>=1){
      return(list(main = cbind(statistic = distribution[1,],pvalue = p_corrected),
                  alpha = alpha))
  }}
  out = list(main = cbind(statistic = distribution[1,],pvalue = p_corrected),
              alpha = alpha)
  return(out)
}

  ############################################################################################################
compute_benjaminin_hochberg <- function(statistic = NULL, pvalue){
  n <- length(pvalue)
  i <- n:1L
  o <- order(pvalue, decreasing = TRUE)
  ro <- order(o)
  pvalue <- pmin(1, cummin(n/i * pvalue[o]))[ro]
  out = list(main = cbind(statistic = statistic,pvalue = pvalue))
  return(out)
}

############################################################################################################


switch_multcomp = function(multcomp,distribution, threshold,aggr_FUN,laterality,E,H,ndh,pvalue,alpha){
  out= list()
  if("maris_oostenveld"%in%multcomp){
    out$maris_oostenveld <- compute_maris_oostenveld(distribution = distribution, threshold = threshold,
                                                                          aggr_FUN = aggr_FUN, laterality = laterality)}
  if("tfce"%in%multcomp){
    out$tfce <- compute_tfce(distribution = distribution, E = E, H = H, ndh = ndh)}
  if("bonferroni"%in%multcomp){
    out$bonferroni <- compute_bonferroni(pvalue = pvalue, statistic = distribution[1,])}
  if("holm"%in%multcomp){
    out$holm <- compute_holm(pvalue = pvalue, statistic = distribution[1,])}
  if("troendle"%in%multcomp){
    out$troendle <- compute_troendle(distribution = distribution, pvalue = pvalue, alpha = alpha, laterality = laterality)}
  if("benjaminin_hochberg"%in%multcomp){
    out$benjaminin_hochberg <- compute_benjaminin_hochberg(pvalue = pvalue, statistic = distribution[1,])}
  return(out)}

#######################################for multcomp output
cluster_table = function(x,...){
  ct = lapply(x, function(effect){
    unique_cluster = unique(effect$maris_oostenveld$main[,3])
    unique_cluster = unique_cluster[unique_cluster!=0]
    if(length(unique_cluster)==0){
      table = "no cluster above the threshold"
      return(table)}
    tab = t(sapply(unique_cluster,function(i){
      cl_select = effect$maris_oostenveld$main[,3] == i
      timepoint = c(1:length(cl_select))[cl_select]
      c(timepoint[1],timepoint[length(timepoint)],
        effect$maris_oostenveld$main[timepoint[1],1],
         effect$maris_oostenveld$main[timepoint[1],2])

    }))
    tab = data.frame(tab)
    colnames(tab) = c("start","end", "cluster mass", "P(>mass)")
    rownames(tab) = unique_cluster
    tab
  })
  class(ct) = "cluster_table"
  ct
}

###### distribution to p_scale

distribution_to_pscale <- function(distribution, test, lateraltiy){
  if(test == "fisher"){
    out <- apply(distribution,2,function(col){compute_all_pvalue(col,laterality = "left")})
  }else if (test == "t"){
    switch(lateraltiy,
           "right" = {out <- apply(distribution,2,function(col){compute_all_pvalue(col,laterality = "left")})},
           "left" = {out <- apply(distribution,2,function(col){compute_all_pvalue(col,laterality = "right")})},
           "bilateral" = {out <- apply(abs(distribution),2,function(col){compute_all_pvalue(col,laterality = "left")})})}
  return(out)
}


