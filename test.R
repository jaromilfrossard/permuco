library(permuco)
library(Matrix)

lf= list.files("R")
for(lfi in lf){source(paste0("R/",lfi))}

## Cluster-mass for repeated measures ANOVA
## Warning : np argument must be greater (recommendation: np>=5000)

attentionshifting_design2 = attentionshifting_design
contrasts(attentionshifting_design2$visibility) = contr.sum
attentionshifting_design2$visibility = as.numeric(attentionshifting_design2$visibility)


npe = 400
t0=proc.time()
signi = attentionshifting_signal[,350:370]
fx <- clusterlm(signi ~ visibility*emotion,
                         data = attentionshifting_design2,
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
