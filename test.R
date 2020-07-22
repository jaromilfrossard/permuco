rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(permuco)
H=function(mat){mat%*%MASS:::ginv(t(mat)%*%mat)%*%t(mat)}

set.seed(1)

mydata = expand.grid(particip = paste0("P",stringr::str_pad(1:6,2,"left","0")),
                     W1 = paste0("w",1:3),
                     W2 = paste0("w",1:2), stringsAsFactors = F)%>%
  mutate(count=if_else(particip=="P01"&W1=="w1"&W2=="w1",4,1))%>%
  uncount(count)%>%
  nest(data=-particip)%>%
  mutate(B1 = if_else(particip%in%c("P01","P02","P03"),"b1","b2"),
         Cov1 = rnorm(n()),
         Cov1 = Cov1-mean(Cov1))%>%
  unnest(data)%>%
  mutate(y = rnorm(n()))%>%
  arrange(particip,W1,W2,B1)


mydata_avg = mydata%>%
  group_by(particip,W1,W2,B1,Cov1)%>%
  summarise(y = mean(y))%>%
  ungroup()%>%
  arrange(particip,W1,W2,B1)


f = y~Cov1*B1*W1*W2+Error(particip/(W1*W2))
#f = y~W1+Error(particip/(W1))

m0 = aovperm(f,data = mydata,np=5)

m1 = aovperm(f,data = mydata_avg,np=5)


m2 = aovperm(f,data = mydata,
             method ="Rd_replic_kheradPajouh_renaud",
             np=5)

m3=aovperm(f,data = mydata_avg,
        method ="Rd_replic_kheradPajouh_renaud",
        np=5)


cbind(wrong=m0$distribution[1,],
      avg=m1$distribution[1,],
      test =m2$distribution[1,],
      test_avg =m3$distribution[1,])


mm <- argi$mm
zid <- argi$mm_id
factors = names(attr(mm,"contrasts"))
eff = colnames(argi$link)
wfact_link <- sapply(strsplit(eff,":"),function(effi){sum(!effi%in%factors)==0})
link_fact <- argi$link[,wfact_link,drop=F]
wfact_mm <- (attr(mm,"assign"))%in%c(0,link_fact[2,drop=F]+link_fact[3,drop=F])
ZE = permuco:::khatrirao(zid,mm[,wfact_mm,drop=F])
mm_gr=cbind(mm[,wfact_mm,drop=F],ZE)
qr_gr = qr(mm_gr)
#hii = diag(qr.fitted(qr_gr,diag(length(argi$y))))
duplic = duplicated.matrix(mm_gr)

###method
assign = attr(mm,"assign")
select_x = assign==argi$i
#select_between = assign%in%args$link[2,]
select_within = assign == (argi$link[3,argi$i])


z = permuco:::khatrirao(a = argi$mm_id, b = mm[,select_within,drop=F])
z = qr.resid(qr(mm),z)
qr_z = qr(z)


qr_d = qr(mm[,!select_x,drop=F])
rdx = qr.resid(qr_d,mm[,select_x,drop=F])

qr_rdx = qr(rdx)
qr_rdz = qr(qr.resid(qr_d,z))
###method

pry = qr.fitted(qr_gr,qr.resid(qr_d,argi$y))
pry[!duplic]

qr_dd = qr(mm[!duplic,!select_x,drop=F])
qr.resid(qr_dd,mm[!duplic,select_x,drop=F])

qr_rdx = qr(qr.resid(qr_d,mm[,select_x,drop=F]))
qr_rdz = qr(qr.resid(qr_d,z))

den = sum(qr.fitted(qr_rdz,pry)^2)/qr_rdz$rank
num = sum(qr.fitted(qr_rdx,pry)^2)/qr_rdx$rank
den


qr.fyi - m3$y


qr_gr = qr(mm_gr)
Q_gr = qr.Q(qr_gr)


qr.fitted(qr_gr,diag(length(argi$y)))-H(mm_gr)

-as.numeric(H(mm_gr)%*%argi$y)


H(mm_gr)%*%argi$y

sum(abs(Q_gr-Q_grc)>1e-10)
qr_gr$rank

rowSums(Q_gr-Q_grc)
sum(abs(qr.fitted(qr_gr,argi$y)-as.numeric(H(mm_gr)%*%argi$y))> 1e-10)

H(mm_gr)- t(Q_gr)%*%Q_gr

diag(H(mm_gr))

qr.Q(qr_gr)%*%t(qr.Q(qr_gr))-H(mm_gr)


DXZE <- cbind(mm,zid,permuco:::khatrirao(mm,zid))

image(Matrix(round(H(mm_gr),6)))

image(Matrix(H(mm)))

rowSums(qr.Q(qr(DXZE))^2)




# mydata = expand.grid(particip = paste0("P",stringr::str_pad(1:6,2,"left","0")),
#                      W1 = paste0("w",1:3),
#                      W2 = paste0("w",1:2), stringsAsFactors = F)%>%
#   mutate(count=if_else(particip=="P01"&W1=="w1"&W2=="w1",4,1))%>%
#   uncount(count)%>%
#   nest(data=-particip)%>%
#   mutate(B1 = if_else(particip%in%c("P01","P02","P03"),"b1","b2"),
#          Cov1 = rnorm(n()),
#          Cov1 = Cov1-mean(Cov1))%>%
#   unnest(data)%>%
#   mutate(y = rnorm(n()))%>%
#   arrange(particip,W1,W2,B1)
#
#
# mydata_avg = mydata%>%
#   group_by(particip,W1,W2,B1,Cov1)%>%
#   summarise(y = mean(y))%>%
#   ungroup()%>%
#   arrange(particip,W1,W2,B1)










# library(permuco)
# library(Matrix)
rm(list=ls())
set.seed(1)
y = rnorm(10)
x = matrix(rnorm(20),ncol=2)


lf= list.files("R")
for(lfi in lf){source(paste0("R/",lfi))}
pm = Pmat(np = 5000, n = 19,type="s",counting = "all")
dim(pm)


n = 15
mat = t(as.matrix(expand.grid(as.data.frame(t(cbind(rep(1,n),rep(-1,n)))))))
rownames(mat)<-NULL
dim(mat)
expand.grid(A = c(1,-1),B = c(1,-1))

ys = permuco::attentionshifting_signal[,(1:50)*16]
df = permuco::attentionshifting_design


lf= list.files("R/")
for(lfi in lf){source(paste0("R/",lfi))}
pm = Pmat(np = 2000, n = nrow(ys),type="s")
plot(hj_p,effect = "visibility")
plot(cl_hj)

cl_hj<-clusterlm(ys~visibility+emotion,df,P=pm,method ="huh_jhun")


hj_p <- aovperm(ys[,1]~visibility+emotion,df,P=pm,method ="huh_jhun")

hj_s <- aovperm(ys[,1]~visibility+emotion,df,P=hj_s$P,method ="huh_jhun")

aovperm(ys[,1]~visibility+emotion,df,type="sign",method ="huh_jhun")

lmperm(ys[,1]~visibility+emotion,df,np=2000,method ="huh_jhun")

aovperm(ys[,1]~visibility*emotion+Error(id/(visibility*emotion)),df,P=pm)
aovperm(ys[,1]~visibility*emotion+Error(id/(visibility*emotion)),df,type="sign")
aovperm(ys[,1]~visibility*emotion+Error(id/(visibility*emotion)),df,P=pm,method ="Rde")
aovperm(ys[,1]~visibility*emotion+Error(id/(visibility*emotion)),df,type="sign",method ="Rde")



param= expand.grid(method = c("freedman_lane","kennedy","huh_jhun","manly",
                              "terBraak","draper_stoneman","dekker"),stringsAsFactors = F)
model_permuco = list()
for(i in 1:nrow(param)){
  set.seed(1)
  print(i)
  model_permuco[[i]]<- clusterlm(ys~visibility+emotion,df,type="sign",method = param$method[i])}

pm = Pmat(np = 2000, n = nrow(ys),type = "permuta")
a = clusterlm(ys~visibility*emotion,df,P=pm)
a$multiple_comparison$visibility$uncorrected$test_info

sa = summary(a)

class(sa[[1]])

attributes(sa[[1]])

a = clusterlm(ys~visibility*emotion,df,type="sign")

attributes(sa[[1]])

pm = Pmat(np = 2000, n = nrow(ys),type = "permuta")
clusterlm(ys~visibility*emotion+Error(id/(visibility*emotion)),df,P=pm,method = "Rd_")
clusterlm(ys~visibility*emotion+Error(id/(visibility*emotion)),df,type="sign",method = "Rd_")
clusterlm(ys~visibility*emotion+Error(id/(visibility*emotion)),df,P=pm,method = "Rde_")
pm = Pmat(np = 2000, n = nrow(ys),type = "signflip")
clusterlm(ys~visibility*emotion+Error(id/(visibility*emotion)),df,P=pm,method = "Rd_")
clusterlm(ys~visibility*emotion+Error(id/(visibility*emotion)),df,P=pm,method = "Rde_")

model_list = list()
param= expand.grid(method = c("freedman_lane","kennedy","huh_jhun","manly",
                              "terBraak","draper_stoneman","dekker"),
                   type=c("permutation","signflip"),stringsAsFactors = F)
for(i in 1:nrow(param)){
  set.seed(1)
  print(i)
  pm = Pmat(np = 2000, n = nrow(ys),type = param$type[i])
  model_list[[i]]<- clusterlm(ys~visibility+emotion,df,P=pm,method = param$method[i])}

i=1
model_list[[8]]
summary(model_permuco[[i]])[[1]]
summary(model_list[[i]])[[1]]
summary(model_list[[i+7]])[[1]]

model_permuco = list()
param = expand.grid(method = c("freedman_lane","kennedy","huh_jhun","manly",
                              "terBraak","draper_stoneman","dekker"),stringsAsFactors = F)
for(i in 1:nrow(param)){
  set.seed(1)
  pm = permuco:::Pmat(np = 2000, n = 10)
  model_permuco[[i]] = permuco:::aovperm(y~x[,1]+x[,2], P = pm, method = param$method[i])}




lf= list.files("R")
for(lfi in lf){source(paste0("R/",lfi))}



model_list = list()
param= expand.grid(method = c("freedman_lane","kennedy","huh_jhun","manly",
            "terBraak","draper_stoneman","dekker"),
            type=c("permutation","signflip"),stringsAsFactors = F)
for(i in 1:nrow(param)){
  set.seed(1)
  pm = Pmat(np = 2000, n = 10,type = param$type[i])
  model_list[[i]] = aovperm(y~x[,1]+x[,2], P = pm, method = param$method[i])}

i =2
model_permuco[[i]]
model_list[[i]]


Pmat_product(x,P = as.matrix(pm[,2]),type = "signflip")
x
pis = lapply(1:np(pm),function(pi){pm[,pi]})


lmp = lmperm(y~x,P = pm)
lmp$table

n = 10
npp =20
pm = Pmat(np = npp, n = n,type = "sign")
Pmat_product(1:n,pm)
pm = Pmat(np = npp, n = n,type = "perm")
Pmat_product(1:8,pm)
pm = Pmat(10,10,type = "signflip","all")
pm[,2,drop=F]

as.matrix(pm)


pm[,2,drop=F]
as.matrix(pm)
attributes(a)



attentionshifting_design2 = attentionshifting_design
contrasts(attentionshifting_design2$visibility) = contr.sum
attentionshifting_design2$visibility = as.numeric(attentionshifting_design2$visibility)

cluster_fisher_myfunp = permuco:::cluster_kennedy
npe = 400
t0=proc.time()
signi = attentionshifting_signal[,350:370]
fx <- clusterlm(signi ~ visibility*emotion,
                         data = attentionshifting_design2,method = "myfunp",new_method=T,
                         np = npe,test="fisher",multcomp = c("clustermass", "troendl","tfce","minP"),return_distribution = T)




plot(fx,multcomp = "troendle",distinctDVs=T)

a = c("12345")
b = c("12")
nchar(b)<-5



d = fx$multiple_comparison$visibility$uncorrected$distribution

compute_minP(d,alternative = "two.sided")

summary(fx,multcomp = "minP",table_type="full")[[1]]
summary(fx,multcomp = "troendle",table_type="full")[[3]]
fx$multiple_comparison$visibility$

rnd <- clusterlm(signi ~ visibility*emotion+Error(id/(visibility*emotion)),
                         data = attentionshifting_design2,
                         np = npe,test="fisher",multcomp = c("clustermass", "troendl","tfce"))
a = summary(fx)
attributes(a$visibility)
electrod_O1$multiple_comparison$visibility$uncorrected$test_info

a = summary(electrod_O1,multcomp = "troendle",table_type = "full")
summary(electrod_O1,multcomp = "troendle")
summary(electrod_O1,table_type = "full",multcomp = "tfce")
full_table(electrod_O1$multiple_comparison)

class(a$visibility)

electrod_O1$multiple_comparison$visibility$uncorrected$test_info
cm =  cluster_table(electrod_O1$multiple_comparison,"troendle")
x = electrod_O1$multiple_comparison

multcomp= "troendle"

cluster_table(electrod_O1$multiple_comparison,"clustermass")
attributes(cm$visibility)
cm$visibility
multcomp="troendle"



proc.time()-t0


t0=proc.time()
electrod_O1p <- clusterlm(attentionshifting_signal ~ visibility*emotion*direction, data = attentionshifting_design2,
                         test="t",multcomp = "troendle2",P=electrod_O1$P)
proc.time()-t0

unique(summary(electrod_O1,multcomp="troendle")[,2])
unique(summary(electrod_O1p,multcomp="troendle2")[,2])

cbind(electrod_O1p$multiple_comparison$visibility$troendle$main[,2],electrod_O1$multiple_comparison$visibility$troendle$main[,2])
sum(summary(electrod_O1,multcomp="troendle")[,2]<0.05)

plot(electrod_O1,multcomp = "troendle")
plot(electrod_O1,alternative = "less",multcomp = "troendle")
plot(electrod_O1,alternative = "greater",multcomp = "troendle")


plot(electrod_O1, nbbaselinepts = 200, nbptsperunit = 1024)
