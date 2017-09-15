cluster_Rd_kheradPajouh_renaud_rnd =function(args){
  #select x
  mm = args$mm
  assign = attr(mm,"assign")
  select_x = assign==args$i
  select_between = assign%in%args$link[2,]
  select_within = assign == (args$link[3,args$i])

  ##zmat
  z = khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
  z = qr.resid(qr(mm),z)
  qr_z = qr(z)

  qr_d = qr(mm[,!select_x,drop=F])
  rdx = qr.resid(qr_d,mm[,select_x,drop=F])


  qr_rdx = qr(rdx)
  qr_rdz = qr(qr.resid(qr_d,z))
  ry = qr.resid(qr_d,args$y)

  #####permutation
  out = t(apply(args$P,2,function(pi){
  colSums(qr.fitted(qr_rdx,ry[pi,,drop=F])^2)/colSums(qr.fitted(qr_rdz,ry[pi,,drop=F])^2)}))*(qr_rdz$rank/qr_rdx$rank)
  return(out)}

##################
cluster_Rde_kheradPajouh_renaud_rnd =function(args){
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
  out = t(apply(args$P,2,function(pi){
    colSums(qr.fitted(qr_rdex,ry[pi,,drop=F])^2)/colSums(qr.fitted(qr_rdez,ry[pi,,drop=F])^2)}))*(qr_rdez$rank/qr_rdex$rank)
  return(out)}

