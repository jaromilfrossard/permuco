###freedman_lane=============================================================
t_freedman_lane <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  #prdy <- matrix(qr.resid(qr_d, args$y)[args$P],ncol = np(args$P))
  prdy <- Pmat_product(qr.resid(qr_d, args$y),args$P)

  #statistic
  qr.coef(qr(rdx), prdy)/sqrt(colSums(qr.resid(qr_mm, prdy)^2)/sum(rdx^2)) * sqrt(length(args$y)-qr_mm$rank)
}

fisher_freedman_lane <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <-qr(rdx)
  prdy <- Pmat_product(qr.resid(qr_d, args$y),args$P)

  #statistic
  colSums(qr.fitted(qr_rdx, prdy)^2)/colSums(qr.resid(qr_mm, prdy)^2)* (length(args$y)-qr_mm$rank)/(qr_rdx$rank)
}



###manly=============================================================
t_manly <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  #py <- matrix(args$y[args$P],ncol = np(args$P))
  py <- Pmat_product(args$y-mean(args$y),args$P)+mean(args$y)

  #statistic
  qr.coef(qr(rdx), py)/sqrt(colSums(qr.resid(qr_mm, py)^2)/sum(rdx^2)) * sqrt(length(args$y)-qr_mm$rank)
}

fisher_manly <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[,select_x, drop = F])
  qr_rdx = qr(rdx)
  #py <- matrix(args$y[args$P],ncol = np(args$P))
  py <- Pmat_product(args$y-mean(args$y),args$P)+mean(args$y)

  #statistic
  colSums(qr.fitted(qr_rdx, py)^2)/colSums(qr.resid(qr_mm, py)^2)* (length(args$y)-qr_mm$rank)/(qr_rdx$rank)
}




###kennedy=============================================================
t_kennedy <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <-qr(rdx)
  #prdy <- matrix(qr.resid(qr_d, args$y)[args$P],ncol = np(args$P))
  prdy <- Pmat_product(qr.resid(qr_d, args$y),args$P)

  #statistic
  qr.coef(qr_rdx, prdy)/sqrt(colSums(qr.resid(qr_rdx, prdy)^2)/sum(rdx^2)) * sqrt(length(args$y)-qr_mm$rank)
}

fisher_kennedy <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <-qr(rdx)
  #prdy <- matrix(qr.resid(qr_d, args$y)[args$P],ncol = np(args$P))
  prdy <- Pmat_product(qr.resid(qr_d, args$y),args$P)

  #statistic
  colSums(qr.fitted(qr_rdx, prdy)^2)/colSums(qr.resid(qr_rdx, prdy)^2)* (length(args$y)-qr_mm$rank)/(qr_rdx$rank)
}




###dekker=============================================================
t_dekker <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  ry <- qr.resid(qr_d, args$y)

  #statistic
  type = attr(args$P,"type")
  apply(args$P,2,function(pi){
    #rprdx = qr.resid(qr_d,rdx[pi,,drop = F])
    rprdx = qr.resid(qr_d, Pmat_product(x = rdx,P = pi,type = type))
    qr_rdprx = qr(rprdx)
    qr.coef(qr_rdprx, ry)[1]/sqrt(sum(qr.resid(qr_rdprx, ry)^2)/sum(rprdx^2))}) *sqrt(length(ry)-qr_mm$rank)
}

fisher_dekker<- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <-qr(rdx)
  ry <- qr.resid(qr_d, args$y)

  #statistic
  type = attr(args$P,"type")
  apply(args$P,2,function(pi){
    #prdx = rdx[pi,,drop = F]
    prdx = Pmat_product(x = rdx,P = pi,type = type)
    qr_rdprx = qr(qr.resid(qr_d,prdx))
    sum(qr.fitted(qr_rdprx, ry)^2)/sum(qr.resid(qr_rdprx, ry)^2)}) * (length(ry)-qr_mm$rank)/(qr_rdx$rank)

  }




###draper stoneman=============================================================
t_draper_stoneman <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])

  #statistic
  type = attr(args$P,"type")
  apply(args$P,2,function(pi){
    #px = args$mm[pi,select_x, drop = F]
    px = Pmat_product(x = args$mm[,select_x,drop = F], P = pi,type = type)

    qr_dpx = qr(cbind(px,args$mm[,!select_x, drop = F]))
    qr.coef(qr_dpx, args$y)[1]/sqrt(sum(qr.resid(qr_dpx, args$y)^2)/sum(qr.resid(qr_d,px)^2))

  }) *sqrt(length(args$y)-qr_mm$rank)
  #qr.coef(qr(rdx), prdy)/sqrt(colSums(qr.resid(qr_mm, prdy)^2)/sum(rdx^2)) * sqrt(length(args$y)-qr_mm$rank)
}

fisher_draper_stoneman<- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  qr_rdx <- qr(qr.resid(qr_d, args$mm[, select_x, drop = F]))

  #statistic
  type = attr(args$P,"type")

  apply(args$P,2,function(pi){
    #px = args$mm[pi,select_x,drop = F]
    px = Pmat_product(x = args$mm[,select_x,drop = F], P = pi,type = type)

    qr_dpx = qr(cbind(px,args$mm[,!select_x, drop = F]))
    sum(qr.fitted(qr(qr.resid(qr_d,px)), args$y)^2)/sum(qr.resid(qr_dpx, args$y)^2)

  }) * (length(args$y)-qr_mm$rank)/(qr_rdx$rank)

}



###ter braack=============================================================
t_terBraak <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])

  #py <- matrix(qr.resid(qr_mm, args$y)[args$P],ncol = np(args$P))+qr.fitted(qr_mm, args$y)
  py <- Pmat_product(qr.resid(qr_mm, args$y),args$P)+ qr.fitted(qr_mm, args$y)


  #statistic
  out = (qr.coef(qr(rdx), py) - qr.coef(qr_mm, args$y)[select_x])/sqrt(colSums(qr.resid(qr_mm, py)^2)/sum(rdx^2)) * sqrt(length(args$y)-qr_mm$rank)
  out[1] = qr.coef(qr_mm, args$y)[select_x]/sqrt(sum(qr.resid(qr_mm, args$y)^2)/sum(rdx^2)) * sqrt(length(args$y)-qr_mm$rank)
  out
}


fisher_terBraak <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx =  qr(rdx)

  #pry <- matrix(qr.resid(qr_mm, args$y)[args$P],ncol = np(args$P))
  pry <- Pmat_product(qr.resid(qr_mm, args$y),args$P)+ qr.fitted(qr_mm, args$y)


  #statistic
  out = colSums(qr.fitted(qr_rdx, pry)^2)/colSums(qr.resid(qr_mm, pry)^2) * (length(args$y)-qr_mm$rank)/(qr_rdx$rank)
  out[1] = sum(qr.fitted(qr_rdx, args$y)^2)/sum(qr.resid(qr_mm, args$y)^2) * (length(args$y)-qr_mm$rank)/(qr_rdx$rank)
  out
}




###huh_jhun=============================================================
t_huh_jhun <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  ###creat random roation from space
  qr_o= qr(args$rnd_rotation[1:(length(args$y)-qr_d$rank),1:(length(args$y)-qr_d$rank)])
  omega = qr.Q(qr_o)%*%diag(sign(diag(qr.R(qr_o))))
  ####create orthogonal subspace
  qcd = qr.Q(qr_d,complete = T)[,-c(1:qr_d$rank),drop=F]
  v = omega%*%t(qcd)

  ###reducing data
  vx <- v%*%(args$mm[,select_x, drop = F])
  qr_vx <-qr(vx)

  #pvy <- matrix((v%*%args$y)[args$P],ncol = np(args$P))
  pvy <- Pmat_product((v%*%args$y),args$P)

  #statistic
  qr.coef(qr_vx, pvy)/sqrt(colSums(qr.resid(qr_vx, pvy)^2)/sum(rdx^2)) * sqrt(length(args$y)-qr_mm$rank)
}

fisher_huh_jhun <- function(args){
  select_x <- c(1:length(attr(args$mm,"assign")))%in%args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  ###creat random roation from space
  qr_o= qr(args$rnd_rotation[1:(length(args$y)-qr_d$rank),1:(length(args$y)-qr_d$rank)])
  omega = qr.Q(qr_o)%*%diag(sign(diag(qr.R(qr_o))))
  ####create orthogonal subspace
  qcd = qr.Q(qr_d,complete = T)[,-c(1:qr_d$rank),drop=F]
  v = omega%*%t(qcd)

  ###reducing data
  vx <- v%*%(args$mm[,select_x, drop = F])
  qr_vx <-qr(vx)

  #pvy <- matrix((v%*%args$y)[args$P],ncol = np(args$P))
  pvy <- Pmat_product((v%*%args$y),args$P)

  #statistic
  colSums(qr.fitted(qr_vx, pvy)^2)/colSums(qr.resid(qr_vx, pvy)^2)* (length(args$y)-qr_mm$rank)/(qr_vx$rank)

  }


