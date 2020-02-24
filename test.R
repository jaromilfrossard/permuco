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
allPerms::

ys = permuco::attentionshifting_signal[,(1:50)*16]
df = permuco::attentionshifting_design


lf= list.files("R")
for(lfi in lf){source(paste0("R/",lfi))}
pm = Pmat(np = 2000, n = nrow(ys),type="s")

aovperm(ys[,1]~visibility+emotion,df,P=pm,method ="huh_jhun")
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
