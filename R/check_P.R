check_P <- function(P, method, test, n, ncol_x2, ncol_x) {
  check_format <- function(p){
    out <- F
    switch(class(p),
           "matrix" = {out <- !is.logical(tryCatch(as.Pmat(P),error = function(e){T}))},
           "Pmat" = {out <- T})
    return(out)}

  nrow_n <- method != "huh_jhun"
  check <- F
  switch(paste(nrow_n),
         "TRUE" = {if(check_format(P)& NROW(P) == n){check <- T}},  ###NOT huh jhun control class et dimension
         "FALSE" = { ####huh jhhun
           switch(test,
                  "t" = {if(check_format(P)& NROW(P) == n-ncol_x+1){check <- T}},
                  "fisher" = {
                    switch(class(P),
                           "list" = {
                             if(all(c(sapply(P,function(p){check_format(p)}),
                                      sapply(P,function(p){NROW(p)})== n - ncol_x + ncol_x2))){check <- T}})})})
  return(check)
}