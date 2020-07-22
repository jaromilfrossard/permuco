



############################################################################################################








############################################################################################################



############################################################################################################



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





########################################################################################




############################################################################################################




#######################################for multcomp output






###### distribution to p_scale



#####

# computed degree of freedom for the fixed effect model
#
# used for the threshold based on the 95 quantile


# computed degree of freedom for the random effect model
#
# used for the threshold based on the 95 quantile



