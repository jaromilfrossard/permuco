# compute_tfce_statistic <- function(signal, E, H, dh, dhi){
#   ##reduce dhi to signal max
#   dhi <- dhi[dhi<=max(signal)]
#   ##compute cluster by dh
#   cl_log <- t(sapply(dhi,function(hi){signal > hi}))
#   cl_beg <- (cl_log-cbind(F,cl_log[,-NCOL(cl_log), drop=F]))==1
#   cl_num <- t(apply(cl_beg,1,cumsum));cl_num[!cl_log]<-0 ##unique integer by cluster by row
#   extend_by_h <- t(apply(cl_num,1,function(cl_row){
#     c(0,sapply(1:max(cl_row),function(i){
#       (sum(cl_row==i))^E
#     }))[cl_row+1]
#   }))
#   ###multiple extend by height
#   statistic <- colSums(extend_by_h*dhi^H)*dh
#   statistic
# }