library(permuco)



## Cluster-mass for repeated measures ANOVA
## Warning : np argument must be greater (recommendation: np>=5000)

attentionshifting_design2 = attentionshifting_design
contrasts(attentionshifting_design2$visibility) = contr.sum
attentionshifting_design2$visibility = as.numeric(attentionshifting_design2$visibility)

electrod_O1 <- clusterlm(attentionshifting_signal ~ visibility*emotion*direction, data = attentionshifting_design2,
                         np = 50,test="t")

electrod_O1$cluster_table_less


plot(electrod_O1, nbbaselinepts = 200, nbptsperunit = 1024)
