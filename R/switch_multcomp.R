switch_multcomp = function(multcomp,distribution, threshold,aggr_FUN,alternative,E,H,ndh,pvalue,alpha,border = border,
                           depth_scale = depth_scale){
  out <- list()
  if("clustermass"%in%multcomp){
    out$clustermass <- compute_clustermass(distribution = distribution, threshold = threshold,
                                           aggr_FUN = aggr_FUN, alternative = alternative)}
  if("clusterdepth"%in%multcomp){
    out$clusterdepth <- compute_clusterdepth(distribution = distribution, threshold = threshold,
                                             alternative = alternative, depth_scale = depth_scale,
                                             border = border)}
  if("clusterdepth_head"%in%multcomp){
    out$clusterdepth <- compute_clusterdepth_head(distribution = distribution, threshold = threshold,
                                             alternative = alternative, depth_scale = depth_scale,
                                             border = border)}
  if("tfce"%in%multcomp){
    out$tfce <- compute_tfce(distribution = distribution, alternative = alternative, E = E, H = H, ndh = ndh)}
  if("bonferroni"%in%multcomp){
    out$bonferroni <- compute_bonferroni(pvalue = pvalue, statistic = distribution[1,])}
  if("holm"%in%multcomp){
    out$holm <- compute_holm(pvalue = pvalue, statistic = distribution[1,])}
  if("troendle"%in%multcomp){
    out$troendle <- compute_troendle(distribution = distribution, alternative = alternative)}
  if("minP"%in%multcomp){
    out$minP <- compute_minP(distribution = distribution, alternative = alternative)}
  if("maxT"%in%multcomp){
    out$maxT <- compute_maxT(distribution = distribution, alternative = alternative)}
  if("stepdownmaxT"%in%multcomp){
    out$stepdownmaxT <- compute_stepdownmaxT(distribution = distribution, alternative = alternative)}
  if("benjamini_hochberg"%in%multcomp){
    out$benjamini_hochberg <- compute_benjamini_hochberg(pvalue = pvalue, statistic = distribution[1,])}
  return(out)}