library(permuco)

lf= list.files("R")
for(lfi in lf){source(paste0("R/",lfi))}

## Cluster-mass for repeated measures ANOVA
## Warning : np argument must be greater (recommendation: np>=5000)

attentionshifting_design2 = attentionshifting_design
contrasts(attentionshifting_design2$visibility) = contr.sum
attentionshifting_design2$visibility = as.numeric(attentionshifting_design2$visibility)
t0=proc.time()

signi = attentionshifting_signal[,400:410]
electrod_O1 <- clusterlm(signi ~ visibility*emotion, data = attentionshifting_design2,
                         np = 10,test="fisher",multcomp = c("clustermass", "troendle"))

a = summary(electrod_O1)
summary(electrod_O1,multcomp = "troendle")

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
