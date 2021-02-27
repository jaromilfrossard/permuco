context("testing permutation methods without random effects")

library(permuco)

set.seed(42)

all_methods_t <-list(
  lm_t_k_p =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="kennedy"),
  lm_t_fl_p =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="freedman_lane"),
  lm_t_d_p =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="dekker"),
  lm_t_hj_p =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="huh_jhun"),
  lm_t_m_p =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="manly"),
  lm_t_tb_p =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="terBraak"),

  lm_t_k_sf =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="kennedy",type ="signflip"),
  lm_t_fl_sf =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="freedman_lane",type ="signflip"),
  lm_t_d_sf =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="dekker",type ="signflip"),
  lm_t_hj_sf =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="huh_jhun",type ="signflip"),
  lm_t_m_sf =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="manly",type ="signflip"),
  lm_t_tb_sf =  lmperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="terBraak",type ="signflip"))

all_methods_f <-list(
  lm_f_k_p =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="kennedy"),
  lm_f_fl_p =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="freedman_lane"),
  lm_f_d_p =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="dekker"),
  lm_f_hj_p =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="huh_jhun"),
  lm_f_m_p =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="manly"),
  lm_f_tb_p =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="terBraak"),

  lm_f_k_sf =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="kennedy",type ="signflip"),
  lm_f_fl_sf =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="freedman_lane",type ="signflip"),
  lm_f_d_sf =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="dekker",type ="signflip"),
  lm_f_hj_sf =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="huh_jhun",type ="signflip"),
  lm_f_m_sf =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="manly",type ="signflip"),
  lm_f_tb_sf =  aovperm(cost ~ sex*insurance, data = emergencycost, np = 20, method ="terBraak",type ="signflip"))


all_methods_t$lm_t_k_p$distribution[1,]

test_that("all methods should give the same t statistics",{
  tstat_distr <- t(sapply(all_methods_t,function(mi){mi$distribution[1,c("sex1","insurance1", "sex1:insurance1")]}))
  tstat_table <- t(sapply(all_methods_t,function(mi){
    si <- summary(mi)
    si[rownames(si)%in%c("sex1","insurance1", "sex1:insurance1"),3]
    }))

  tstat <- round(rbind(tstat_distr,tstat_table),12)
  expect_setequal(tstat[,1], tstat[1,1])
  expect_setequal(tstat[,2], tstat[1,2])
  expect_setequal(tstat[,3], tstat[1,3])
})



test_that("all methods should give the same F statistics",{
  fstat_distr <- t(sapply(all_methods_f,function(mi){mi$distribution[1,c("sex","insurance", "sex:insurance")]}))
  fstat_table <- t(sapply(all_methods_f,function(mi){
    si <- summary(mi)
    si[rownames(si)%in%c("sex","insurance", "sex:insurance"),3]
  }))

  fstat <- round(rbind(fstat_distr,fstat_table),12)


  fstat <- round(fstat,12)
  expect_setequal(fstat[,1], fstat[1,1])
  expect_setequal(fstat[,2], fstat[1,2])
  expect_setequal(fstat[,3], fstat[1,3])
})


