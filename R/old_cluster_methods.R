








# #' @export
# summary.cluster_table <- function(object,...){
#   object
# }


# summary_multcomp <- function(object, multcomp, alternative){
#   if(!(multcomp %in% object$multcomp)){
#     stop(paste("The choosen multiple comparisons procedure does not match with the ones computed in the object.
#                   Choose between ", paste(object$multcomp,sep=" "),".", sep=""))
#   }
#   switch(alternative,
#          "two.sided" = {mc = object$multiple_comparison},
#          "greater" = {mc = object$multiple_comparison_greater},
#          "less" = {mc = object$multiple_comparison_less})
#
#   out = lapply(seq_along(mc),function(i){
#     out = mc[[i]][[multcomp]]$main[,1:2,drop = F]
#     colnames(out) = paste(names(mc)[i], colnames(out))
#     out})
#   out = do.call("cbind",out)
#   return(out)
# }







