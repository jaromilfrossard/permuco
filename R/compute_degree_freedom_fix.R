compute_degree_freedom_fix = function(test,mm,assigni){
  qr_mm <- qr(mm)
  switch(test,
         "t" = {rep(NROW(mm)-qr_mm$rank,length(assigni))},
         "fisher" = {cbind(dfn = as.numeric(table(assigni[assigni!=0])), dfd =NROW(mm)-qr_mm$rank)})
}