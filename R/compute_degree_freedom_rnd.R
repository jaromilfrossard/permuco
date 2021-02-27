compute_degree_freedom_rnd = function(test = "fisher",mm,mm_id,assigni,link){
  qr_mm <- qr(mm)
  switch(test,
         "t" = {},
         "fisher" = {
           cbind(dfn <- as.numeric(table(assigni[assigni!=0])),
                 dfd <- sapply(1:max(assigni[assigni!=0]),function(i){
                   qr(qr.resid(qr_mm,khatrirao(mm_id, mm[,attr(mm,"assign")==link[3,i],drop=F])))$rank}))
         })
}