context("testing permutation methods without random effects")

library(permuco)

set.seed(42)


all_methods_rnd <-list(
  aov_rd_p = aovperm(iapa ~ bmi*time+ Error(id/(time)),
                      data = jpah2016, method = "Rd_kheradPajouh_renaud",
                      np = 20, type="permutation"),
  aov_rd_sf = aovperm(iapa ~ bmi*time+ Error(id/(time)),
                      data = jpah2016, method = "Rd_kheradPajouh_renaud",
                      np = 20, type="signflip"),
  aov_rde_p = aovperm(iapa ~ bmi*time+ Error(id/(time)),
                      data = jpah2016, method = "Rde_kheradPajouh_renaud",
                      np = 20, type="permutation"),
  aov_rde_sf = aovperm(iapa ~ bmi*time+ Error(id/(time)),
                      data = jpah2016, method = "Rde_kheradPajouh_renaud",
                      np = 20, type="signflip"))

test_that("all methods should give the same F statistics in ranova",{
  fstat_rnd_distr <- t(sapply(all_methods_rnd,function(mi){mi$distribution[1,c("bmi","time","bmi:time")]}))
  fstat_rnd_table <- t(sapply(all_methods_rnd,function(mi){
    si <- summary(mi)
    si[rownames(si)%in%c("bmi","time","bmi:time"),7]
  }))

  fstat_rnd <- round(rbind(fstat_rnd_distr,fstat_rnd_table),12)


  fstat_rnd <- round(fstat_rnd,12)
  expect_setequal(fstat_rnd[,1], fstat_rnd[1,1])
  expect_setequal(fstat_rnd[,2], fstat_rnd[1,2])
  expect_setequal(fstat_rnd[,3], fstat_rnd[1,3])

})