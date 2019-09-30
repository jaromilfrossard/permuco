############################################################################################################
compute_clustermass <-  function(distribution, threshold, aggr_FUN, alternative = "greater"){
    switch(alternative,
           "greater" = {
             threshold <- abs(threshold)
             selected <- distribution > threshold
             extreme = function(x)max(x,na.rm = T)
           },
           "less" = {
             threshold <- -abs(threshold)
             selected <- distribution < threshold
             extreme = function(x)max(x,na.rm = T)},
           "two.sided" = {
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
    compute_pvalue(stat = mi,distribution = mass_distribution, alternative = "greater"))
  pvalue = c(NA,pvalue)[cl[1,]+1]
  out = list(main = cbind(statistic = statistic, pvalue = pvalue, cluster_id = cl[1,]),
              distribution = mass_distribution,
             threshold = threshold)
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
# compute_troendle_old <- function(distribution, pvalue, alpha, alternative){
#   distribution_rank = apply(distribution,2,function(col){compute_all_pvalue(col,alternative = alternative)})
#
#   p_corrected <- rep(1,length(pvalue))
#   rank_uncorr <- rank(pvalue)
#   #loop to recompute minimal value
#   for(urank in sort(unique(rank_uncorr))){
#     which_test <- which(urank==rank_uncorr)
#     pvali <- distribution_rank[,which(urank<=rank_uncorr)]
#     distr_min <- apply(pvali,1,min)
#     p_corri <- compute_pvalue(distribution = distr_min, stat = matrix(pvalue[which_test],nrow=1), alternative = "less")
#     p_corrected[which_test] <- p_corri
#     if(sum(p_corri > alpha)>=1){
#       return(list(main = cbind(statistic = distribution[1,],pvalue = p_corrected),
#                   alpha = alpha))
#   }}
#   out = list(main = cbind(statistic = distribution[1,],pvalue = p_corrected),
#               alpha = alpha)
#   return(out)
# }
compute_troendle = function (distribution, alternative) {
  distribution_rank = apply(distribution,2,function(col){compute_all_pvalue(col,alternative = alternative)})

  pvalue = distribution_rank[1,]

  rank_uncorr <- rank(pvalue)


  # order test from low to high
  test_order = sort(unique(rank_uncorr))
  order_test_list = lapply(1:length(test_order),function(testi){
    ri = test_order[testi]
    list(infos = cbind(order=testi,col=which(rank_uncorr==ri)),
         mins = apply(distribution_rank[,which(rank_uncorr==ri),drop=F],1,min))
  })

  # take cumulative minimum of the ordered distributions of p's

  cummins = t(apply(do.call("cbind",lapply(order_test_list,function(oi)(oi$mins))),
                    1,function(coli){rev(cummin(rev(coli)))}))

  # Compute p-value

  p_corrected = apply(cummins,2,
                      function(coli)compute_pvalue(coli,alternative = "less"))

  # Compute inverse order

  infos = do.call("rbind",lapply(order_test_list,function(oi)(oi$infos)))

  # Inverse order of the increasing p-values.
  p_corrected = (cummax(p_corrected) [infos[,1]])[order(infos[,2])]

  out = list(main = cbind(statistic = distribution[1, ], pvalue = p_corrected))
  return(out)
}



  ############################################################################################################
compute_benjamini_hochberg <- function(statistic = NULL, pvalue){
  n <- length(pvalue)
  i <- n:1L
  o <- order(pvalue, decreasing = TRUE)
  ro <- order(o)
  pvalue <- pmin(1, cummin(n/i * pvalue[o]))[ro]
  out = list(main = cbind(statistic = statistic,pvalue = pvalue))
  return(out)
}

############################################################################################################


switch_multcomp = function(multcomp,distribution, threshold,aggr_FUN,alternative,E,H,ndh,pvalue,alpha){
  out= list()
  if("clustermass"%in%multcomp){
    out$clustermass <- compute_clustermass(distribution = distribution, threshold = threshold,
                                                                          aggr_FUN = aggr_FUN, alternative = alternative)}
  if("tfce"%in%multcomp){
    out$tfce <- compute_tfce(distribution = distribution, E = E, H = H, ndh = ndh)}
  if("bonferroni"%in%multcomp){
    out$bonferroni <- compute_bonferroni(pvalue = pvalue, statistic = distribution[1,])}
  if("holm"%in%multcomp){
    out$holm <- compute_holm(pvalue = pvalue, statistic = distribution[1,])}
  if("troendle"%in%multcomp){
    out$troendle <- compute_troendle(distribution = distribution, alternative = alternative)}
  if("benjamini_hochberg"%in%multcomp){
    out$benjamini_hochberg <- compute_benjamini_hochberg(pvalue = pvalue, statistic = distribution[1,])}
  return(out)}

#######################################for multcomp output
cluster_table = function(x,...){
  dotargs = list(...)
  ct = lapply(1:length(x), function(j){
    effect = x[[j]]
    unique_cluster = unique(effect$clustermass$main[,3])
    unique_cluster = unique_cluster[unique_cluster!=0]
    if(length(unique_cluster)==0){
      attr(table,"effect_name") = names(x)[j]
      attr(table,"threshold") = effect$clustermass$threshold
      table = paste(names(x)[j], ", no cluster above a threshold of : ", round(effect$clustermass$threshold, 5),
                    sep= "")
      return(table)}
    tab = t(sapply(unique_cluster,function(i){
      cl_select = effect$clustermass$main[,3] == i
      timepoint = c(1:length(cl_select))[cl_select]
      c(timepoint[1],timepoint[length(timepoint)],
        effect$clustermass$main[timepoint[1],1],
         effect$clustermass$main[timepoint[1],2])

    }))
    tab = data.frame(tab)
    colnames(tab) = c("start","end", "cluster mass", "P(>mass)")
    rownames(tab) = unique_cluster
    attr(tab,"threshold") = effect$clustermass$threshold
    attr(tab,"effect_name") = names(x)[j]
    class(tab) = append("cluster_table",class(tab))
    tab
  })
  class(ct) = append("listof_cluster_table",class(ct))
  names(ct) = names(x)
  ct
}

###### distribution to p_scale

distribution_to_pscale <- function(distribution, test, alternative){
  if(test == "fisher"){
    out <- apply(distribution,2,function(col){compute_all_pvalue(col,alternative = "less")})
  }else if (test == "t"){
    switch(alternative,
           "greater" = {out <- apply(distribution,2,function(col){compute_all_pvalue(col,alternative = "less")})},
           "less" = {out <- apply(distribution,2,function(col){compute_all_pvalue(col,alternative = "greater")})},
           "two.sided" = {out <- apply(abs(distribution),2,function(col){compute_all_pvalue(col,alternative = "less")})})}
  return(out)
}

#####

# computed degree of freedom for the fixed effect model
#
# used for the threshold based on the 95 quantile
compute_degree_freedom_fix = function(test,mm,assigni){
  qr_mm = qr(mm)
  switch(test,
         "t" = {rep(NROW(mm)-qr_mm$rank,length(assigni))},
         "fisher" = {cbind(dfn = as.numeric(table(assigni[assigni!=0])), dfd =NROW(mm)-qr_mm$rank)})
}

# computed degree of freedom for the random effect model
#
# used for the threshold based on the 95 quantile
compute_degree_freedom_rnd = function(test = "fisher",mm,mm_id,assigni,link){
  qr_mm = qr(mm)
  switch(test,
         "t" = {},
         "fisher" = {
           cbind(dfn = as.numeric(table(assigni[assigni!=0])),
                 dfd=sapply(1:max(assigni[assigni!=0]),function(i){
                   qr(qr.resid(qr_mm,khatrirao(mm_id, mm[,attr(mm,"assign")==link[3,i],drop=F])))$rank}))
         })
}


