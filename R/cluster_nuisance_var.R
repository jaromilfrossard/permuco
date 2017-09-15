#freedman lane================================
cluster_freedman_lane <- function(args){
  ##test selection
  switch(args$test,
         "fisher"= {funT = function(qr_rdx, qr_mm, prdy){
           colSums(qr.fitted(qr_rdx, prdy)^2)/colSums(qr.resid(qr_mm, prdy)^2)* (NROW(prdy)-qr_mm$rank)/(qr_rdx$rank)}
         },
         "t" = {funT = function(qr_rdx, qr_mm, prdy){
           as.numeric(qr.coef(qr_rdx, prdy))/sqrt(colSums(qr.resid(qr_mm, prdy)^2)/sum(rdx^2)) * sqrt(NROW(args$y)-qr_mm$rank)}
         })

  ##effect selection
  select_x <- c(1:length(attr(args$mm,"assign"))) %in% args$colx

  ##data reduction
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <- qr(rdx)
  rdy <- qr.resid(qr_d, args$y)


  out = apply(args$P,2,function(pi){
    funT(qr_rdx = qr_rdx, qr_mm = qr_mm, prdy = rdy[pi,, drop = F])})
  return(out)}


#kennedy================================
cluster_kennedy <- function(args){
  ##test selection
  switch(args$test,
         "fisher"= {funT = function(qr_rdx, qr_mm, prdy){
           colSums(qr.fitted(qr_rdx, prdy)^2)/colSums(qr.resid(qr_rdx, prdy)^2)* (NROW(prdy)-qr_mm$rank)/(qr_rdx$rank)}
         },
         "t" = {funT = function(qr_rdx, qr_mm, prdy){
           as.numeric(qr.coef(qr_rdx, prdy))/sqrt(colSums(qr.resid(qr_rdx, prdy)^2)/sum(rdx^2)) * sqrt(NROW(args$y)-qr_mm$rank)}
         })

  ##effect selection
  select_x <- c(1:length(attr(args$mm,"assign"))) %in% args$colx

  ##data reduction
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <- qr(rdx)
  rdy <- qr.resid(qr_d, args$y)


  out = apply(args$P,2,function(pi){
    funT(qr_rdx = qr_rdx, qr_mm = qr_mm, prdy = rdy[pi,, drop = F])})
  return(out)}



#terBraak================================
cluster_terBraak <- function(args){
  ##test selection
  switch(args$test,
         "fisher"= {funT = function(qr_rdx, qr_mm, pry){
           colSums(qr.fitted(qr_rdx, pry)^2)/colSums(qr.resid(qr_mm, pry)^2)* (NROW(pry)-qr_mm$rank)/(qr_rdx$rank)}
         },
         "t" = {funT = function(qr_rdx, qr_mm, pry){
           as.numeric(qr.coef(qr_rdx, pry))/sqrt(colSums(qr.resid(qr_mm, pry)^2)/sum(rdx^2)) * sqrt(NROW(args$y)-qr_mm$rank)}
         })



  ##effect selection
  select_x <- c(1:length(attr(args$mm,"assign"))) %in% args$colx

  ##data reduction
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <- qr(rdx)
  rdy <- qr.resid(qr_d, args$y)
  rmmy <- qr.resid(qr_mm, args$y)


  out = apply(args$P,2,function(pi){
    funT(qr_rdx = qr_rdx, qr_mm = qr_mm, pry = rmmy[pi,, drop = F])})

  out[,1] = funT(qr_rdx = qr_rdx, qr_mm = qr_mm, pry = rdy)



  return(out)}





#manly================================
cluster_manly <- function(args){
  ##test selection
  switch(args$test,
         "fisher"= {funT = function(qr_rdx, qr_mm, py){
           colSums(qr.fitted(qr_rdx, py)^2)/colSums(qr.resid(qr_mm, py)^2)* (NROW(py)-qr_mm$rank)/(qr_rdx$rank)}
         },
         "t" = {funT = function(qr_rdx, qr_mm, py){
           as.numeric(qr.coef(qr_rdx, py))/sqrt(colSums(qr.resid(qr_mm, py)^2)/sum(rdx^2)) * sqrt(NROW(args$y)-qr_mm$rank)}
         })

  ##effect selection
  select_x <- c(1:length(attr(args$mm,"assign"))) %in% args$colx

  ##data reduction
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <- qr(rdx)


  out = apply(args$P,2,function(pi){
    funT(qr_rdx = qr_rdx, qr_mm = qr_mm, py = args$y[pi,, drop = F])})
  return(out)}




#draperstoneman================================
cluster_draper_stoneman <- function(args){
  ##test selection
  switch(args$test,
         "fisher"= {funT = function(qr_rdpx, qr_pmm, y, qr_mm, qr_rdx, rdpx){
           colSums(qr.fitted(qr_rdpx, y)^2)/colSums(qr.resid(qr_pmm, y)^2)* (NROW(y)-qr_mm$rank)/(qr_rdx$rank)}
         },
         "t" = {funT = function(qr_rdpx, qr_pmm, y, qr_mm, qr_rdx, rdpx){
           as.numeric(qr.coef(qr_rdpx, y))/sqrt(colSums(qr.resid(qr_pmm, y)^2)/sum(rdpx^2)) * sqrt(NROW(y)-qr_mm$rank)}
         })

  ##effect selection
  select_x <- c(1:length(attr(args$mm,"assign"))) %in% args$colx

  ##data reduction
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <- qr(rdx)

  out = apply(args$P,2,function(pi){
    rdpx = qr.resid(qr_d,args$mm[pi,select_x, drop=F])
    qr_rdpx = qr(rdpx)
    qr_pmm = qr(cbind(args$mm[,!select_x, drop=F],args$mm[pi,select_x, drop=F]))
    funT(qr_rdpx = qr_rdpx, qr_pmm = qr_pmm, y = args$y, qr_mm = qr_mm, qr_rdx = qr_rdx, rdpx = rdpx)})
  return(out)}


#dekker================================

cluster_dekker <- function(args){
  ##test selection
  switch(args$test,
         "fisher"= {funT = function(qr_rdprx, ry, qr_mm,qr_rdx,rdprx){
           colSums(qr.fitted(qr_rdprx, ry)^2)/colSums(qr.resid(qr_rdprx, ry)^2)* (NROW(ry)-qr_mm$rank)/(qr_rdx$rank)}
         },
         "t" = {funT = function(qr_rdprx, ry, qr_mm,qr_rdx,rdprx){
           as.numeric(qr.coef(qr_rdprx, ry))/sqrt(colSums(qr.resid(qr_rdprx, ry)^2)/sum(rdprx^2)) * sqrt(NROW(ry)-qr_mm$rank)}
         })

  ##effect selection
  select_x <- c(1:length(attr(args$mm,"assign"))) %in% args$colx

  ##data reduction
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <- qr(rdx)
  ry = qr.resid(qr_d,args$y)

  out = apply(args$P,2,function(pi){
    rdprx = qr.resid(qr_d,rdx[pi,, drop=F])
    qr_rdprx = qr(rdprx)
    #qr_pmm = qr(cbind(args$mm[,!select_x, drop=F],rdx[pi,, drop=F]))

    funT(qr_rdprx = qr_rdprx, ry = ry, qr_mm = qr_mm, qr_rdx = qr_rdx, rdprx = rdprx)})
  return(out)}




#huh_jhun================================
cluster_huh_jhun <- function(args){
  ##test selection
  switch(args$test,
         "fisher"= {funT = function(qr_vx, qr_mm, pvy, rdx){
           colSums(qr.fitted(qr_vx, pvy)^2)/colSums(qr.resid(qr_vx, pvy)^2)* (NROW(args$y)-qr_mm$rank)/(qr_vx$rank)}
         },
         "t" = {funT = function(qr_vx, qr_mm, pvy, rdx){
           as.numeric(qr.coef(qr_vx, pvy))/sqrt(colSums(qr.resid(qr_vx, pvy)^2)/sum(rdx^2)) * sqrt(NROW(args$y)-qr_mm$rank)}
         })

  ##effect selection
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])

  ###creat random roation from space
  qr_o= qr(args$rnd_rotation[1:(NROW(args$y)-qr_d$rank),1:(NROW(args$y)-qr_d$rank)])
  omega = qr.Q(qr_o)%*%diag(sign(diag(qr.R(qr_o))))

  ####create orthogonal subspace
  qcd = qr.Q(qr_d,complete = T)[,-c(1:qr_d$rank),drop=F]
  v = omega%*%t(qcd)



  ###reducing data
  vx <- v%*%(args$mm[,select_x, drop = F])
  qr_vx <-qr(vx)
  vy <- v%*%args$y



  out = apply(args$P,2,function(pi){
    funT(qr_vx = qr_vx, qr_mm = qr_mm, pvy = vy[pi,,drop = F], rdx = rdx)})


  return(out)}

