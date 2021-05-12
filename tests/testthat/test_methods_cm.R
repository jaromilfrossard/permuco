context("testing permutation methods in clustermass")
library(permuco)


allmc <- c("clustermass","clusterdepth","clusterdepth_head","troendle","tfce","minP","maxT",
           "stepdownmaxT","benjamini_hochberg","bonferroni","holm")

all_methods_f <- list(
  cm_f_k_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                      method ="kennedy", multcomp = allmc),
  cm_f_k_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                      method ="kennedy", multcomp = allmc,type = "signflip"),

  cm_f_fl_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                      method ="freedman_lane", multcomp = allmc),
  cm_f_fl_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                       method ="freedman_lane", multcomp = allmc,type = "signflip"),

  cm_f_hj_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                       method ="huh_jhun", multcomp = allmc),
  cm_f_hj_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                        method ="huh_jhun", multcomp = allmc,type = "signflip"),

  cm_f_d_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                       method ="dekker", multcomp = allmc),
  cm_f_d_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                        method ="dekker", multcomp = allmc,type = "signflip"),

  cm_f_m_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                      method ="manly", multcomp = allmc),
  cm_f_m_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                       method ="manly", multcomp = allmc,type = "signflip"),
  cm_f_tb_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                      method ="terBraak", multcomp = allmc),
  cm_f_tb_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                       method ="terBraak", multcomp = allmc,type = "signflip"))




all_methods_t <- list(
  cm_t_k_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                        method ="kennedy", multcomp = allmc,test="t"),
  cm_t_k_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                         method ="kennedy", multcomp = allmc,type = "signflip",test="t"),

  cm_t_fl_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                         method ="freedman_lane", multcomp = allmc,test="t"),
  cm_t_fl_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                          method ="freedman_lane", multcomp = allmc,type = "signflip",test="t"),

  cm_t_hj_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                         method ="huh_jhun", multcomp = allmc,test="t"),
  cm_t_hj_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                          method ="huh_jhun", multcomp = allmc,type = "signflip",test="t"),

  cm_t_d_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                        method ="dekker", multcomp = allmc,test="t"),
  cm_t_d_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                         method ="dekker", multcomp = allmc,type = "signflip",test="t"),

  cm_t_m_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                        method ="manly", multcomp = allmc,test="t"),
  cm_t_m_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                         method ="manly", multcomp = allmc,type = "signflip",test="t"),
  cm_t_tb_p =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                         method ="terBraak", multcomp = allmc,test="t"),
  cm_t_tb_sf =  clusterlm(attentionshifting_signal[,100:120] ~ direction*STAIS_state, data = attentionshifting_design, np = 20,
                          method ="terBraak", multcomp = allmc,type = "signflip",test="t"))






test_that("all methods of clusterlm should give the same F statistics",{
  all_f <- t(sapply(all_methods_f,function(mi)do.call("rbind",lapply(mi$multiple_comparison,function(ei)ei$uncorrected$main[,1,drop=F]))))
  all_f <- round(all_f,11)
  row.names(all_f )=NULL
  for(i in seq_len(ncol(all_f))){expect_setequal(all_f[,i], all_f[1,i])}
  })


test_that("all methods of clusterlm should give the same t statistics",{
  all_t <- lapply(all_methods_t,function(mi){
    cbind(
      do.call("rbind",lapply(mi$multiple_comparison[c("direction1","STAIS_state","direction1:STAIS_state")],function(ei)ei$uncorrected$main[,1,drop=F])),
      do.call("rbind",lapply(mi$multiple_comparison_greater[c("direction1","STAIS_state","direction1:STAIS_state")],function(ei)ei$uncorrected$main[,1,drop=F])),
      do.call("rbind",lapply(mi$multiple_comparison_less[c("direction1","STAIS_state","direction1:STAIS_state")],function(ei)ei$uncorrected$main[,1,drop=F])))
    })
  all_t <- t(do.call("cbind",all_t))
  all_t <- round(all_t,12)
  row.names(all_t )=NULL
  for(i in seq_len(ncol(all_t))){expect_setequal(all_t[,i], all_t[1,i])}
})


test_that("dimension of cluster table should be equal",{
  dim_f <- t(sapply(all_methods_f,function(mi) as.numeric(sapply(summary(mi),dim))))
  row.names(dim_f )=NULL
  dim_t<- t(sapply(all_methods_t,function(mi) as.numeric(sapply(summary(mi)[c("direction1","STAIS_state","direction1:STAIS_state")],dim))))
  row.names(dim_t )=NULL
  dim_t_less<- t(sapply(all_methods_t,function(mi) as.numeric(sapply(summary(mi, alternative ="less")[c("direction1","STAIS_state","direction1:STAIS_state")],dim))))
  row.names(dim_t_less )=NULL
  dim_t_greater<- t(sapply(all_methods_t,function(mi) as.numeric(sapply(summary(mi, alternative ="greater")[c("direction1","STAIS_state","direction1:STAIS_state")],dim))))
  row.names(dim_t_greater )=NULL

  for(i in seq_len(ncol(dim_f))){expect_setequal(dim_f[,i], dim_f[1,i])}
  for(i in seq_len(ncol(dim_t))){expect_setequal(dim_t[,i], dim_t[1,i])}
  for(i in seq_len(ncol(dim_t_less))){expect_setequal(dim_t_less[,i], dim_t_less[1,i])}
  for(i in seq_len(ncol(dim_t_greater))){expect_setequal(dim_t_greater[,i], dim_t_greater[1,i])}

})


