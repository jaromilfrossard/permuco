context("testing permutation methods in clustermass")
library(permuco)


allmc <- c("clustermass","clusterdepth","troendle","tfce","minP","benjamini_hochberg","bonferroni","holm")



all_methods_f <- list(
  cm_f_rd_p =  clusterlm(attentionshifting_signal[,100:120] ~visibility*emotion+Error(id/(visibility*emotion)), data = attentionshifting_design, np = 20,
                        method ="Rd_kheradPajouh_renaud", multcomp = allmc),
  cm_f_rd_sf =  clusterlm(attentionshifting_signal[,100:120] ~visibility*emotion+Error(id/(visibility*emotion)), data = attentionshifting_design, np = 20,
                         method ="Rd_kheradPajouh_renaud", multcomp = allmc,type = "signflip"),
  cm_f_rde_p =  clusterlm(attentionshifting_signal[,100:120] ~visibility*emotion+Error(id/(visibility*emotion)), data = attentionshifting_design, np = 20,
                         method ="Rde_kheradPajouh_renaud", multcomp = allmc),
  cm_f_rde_sf =  clusterlm(attentionshifting_signal[,100:120] ~visibility*emotion+Error(id/(visibility*emotion)), data = attentionshifting_design, np = 20,
                          method ="Rde_kheradPajouh_renaud", multcomp = allmc,type = "signflip")

  )


test_that("all methods of clusterlm should give the same F statistics of ranova",{
  all_f <- t(sapply(all_methods_f,function(mi)do.call("rbind",lapply(mi$multiple_comparison,function(ei)ei$uncorrected$main[,1,drop=F]))))
  all_f <- round(all_f,11)
  row.names(all_f )=NULL
  for(i in seq_len(ncol(all_f))){expect_setequal(all_f[,i], all_f[1,i])}
})


test_that("dimension of cluster table should be equal in ranova",{
  dim_f <- t(sapply(all_methods_f,function(mi) as.numeric(sapply(summary(mi),dim))))
  row.names(dim_f )=NULL
  for(i in seq_len(ncol(dim_f))){expect_setequal(dim_f[,i], dim_f[1,i])}

  })
