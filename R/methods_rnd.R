
fisher_Rd_kheradPajouh_renaud_rnd =function(args){
  #select x
  mm = args$mm
  assign = attr(mm,"assign")
  select_x = assign==args$i
  #select_between = assign%in%args$link[2,]
  select_within = assign == (args$link[3,args$i])

  ##zmat
  z = khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
  z = qr.resid(qr(mm),z)
  qr_z = qr(z)

  qr_d = qr(mm[,!select_x,drop=F])
  rdx = qr.resid(qr_d,mm[,select_x,drop=F])

  qr_rdx = qr(rdx)
  qr_rdz = qr(qr.resid(qr_d,z))

  #####permutation
  #pry = matrix(qr.resid(qr_d,args$y)[args$P], nrow = length(args$y))
  pry <- Pmat_product(qr.resid(qr_d,args$y),args$P)

  den = colSums(qr.fitted(qr_rdz,pry)^2)/qr_rdz$rank
  num = colSums(qr.fitted(qr_rdx,pry)^2)/qr_rdx$rank

  return(c(num/den))

}


fisher_Rde_kheradPajouh_renaud_rnd =function(args){
  ##matrix
  mm = args$mm
  assign = attr(mm,"assign")
  select_x = assign==args$i
  select_within = assign == (args$link[3,args$i])

  within0 = unique(args$link[3,])
  select_within_e = assign%in%(within0[!(within0%in%args$link[3,args$i])])

  ##
  qr_mm = qr(mm)

  ##zmat
  z = khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
  z = qr.resid(qr_mm ,z)
  qr_z = qr(z)

  ##emat
  ee = khatrirao(a = args$mm_id, b = mm[,select_within_e,drop=F])
  ee = qr.resid(qr_mm,ee)
  ee = qr.resid(qr_z,ee)

  qr_d = qr(mm[,!select_x,drop=F])
  qr_de = qr(cbind(mm[,!select_x,drop=F],ee))

  rdex = qr.resid(qr_de,mm[,select_x,drop=F])

  qr_rdex = qr(rdex)
  qr_rdez = qr(qr.resid(qr_de,z))

  ry = qr.resid(qr_de,args$y)

  #####permutation
  #pry = matrix(ry[args$P], nrow = length(ry))
  pry <- Pmat_product(ry,args$P)

  den = colSums(qr.fitted(qr_rdez,pry)^2)/qr_rdez$rank
  num = colSums(qr.fitted(qr_rdex,pry)^2)/qr_rdex$rank

  return(c(num/den))

}


fisher_Rd_replic_kheradPajouh_renaud_rnd <- function(args){
  #select x
  mm <- args$mm
  assign = attr(mm,"assign")
  select_x = assign==args$i
  #select_between = assign%in%args$link[2,]
  select_within = assign == (args$link[3,args$i])


  #####
  zid <- args$mm_id
  factors = names(attr(mm,"contrasts"))
  eff = colnames(args$link)
  wfact_link <- sapply(strsplit(eff,":"),function(effi){sum(!effi%in%factors)==0})
  link_fact <- args$link[,wfact_link,drop=F]
  wfact_mm <- (attr(mm,"assign"))%in%c(0,link_fact[2,,drop=F]+link_fact[3,,drop=F])
  ZE = khatrirao(zid,mm[,wfact_mm,drop=F])
  mm_gr=cbind(mm[,wfact_mm,drop=F],ZE)
  qr_gr = qr(mm_gr)
  hii = diag(qr.fitted(qr_gr,diag(length(args$y))))
  duplic = duplicated.matrix(mm_gr)

  ######


  ##zmat
  z = khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
  z = qr.resid(qr(mm),z)
  qr_z = qr(z)

  qr_d = qr(mm[,!select_x,drop=F])
  rdx = qr.resid(qr_d,mm[,select_x,drop=F])

  qr_rdx = qr(rdx)
  qr_rdz = qr(qr.resid(qr_d,z))

  #####permutation
  #pry = matrix(qr.resid(qr_d,args$y)[args$P], nrow = length(args$y))
  pry = Pmat_product(qr.resid(qr_d,args$y),args$P)
  fpry <-qr.fitted(qr_gr,pry)

  #pry <-(hii)*pry

  qr_rdx = qr(qr.resid(qr_d,mm[,select_x,drop=F]))
  qr_rdz = qr(qr.resid(qr_d,z))


  ###
  fpry_dup = fpry[!duplic,,drop=F]
  qr_d_dup = qr(mm[!duplic,!select_x,drop=F])
  rdx_dup = qr.resid(qr_d_dup,mm[!duplic,select_x,drop=F])
  qr_rdx_dup = qr(rdx_dup)


  z_dup = khatrirao(a = args$mm_id[!duplic,,drop=F], b = mm[!duplic,select_within,drop=F])
  z_dup  = qr.resid(qr(mm[!duplic,,drop=F]),z_dup )
  qr_rdz_dup = qr(qr.resid(qr_d_dup,z_dup))
  ###


  den = colSums((qr.fitted(qr_rdz_dup,fpry_dup))^2)/qr_rdz_dup$rank
  num = colSums((qr.fitted(qr_rdx_dup,fpry_dup))^2)/qr_rdx_dup$rank

  return(c(num/den))

}





