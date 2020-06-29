require(glmmTMB)
require(tidyr)
require(ggplot2)
require(doBy)
require(DHARMa)
require(broomExtra)
require(dotwhisker)
require(broom)
require(data.table)
require(sjPlot)
require(e1071)
require(insight)
require(dplyr)
require(Hmisc)
require(ggeffects)

#==============================================================================================#
###READ IN DATA
seed1 <- read.csv(file="seedmaster_climNA.csv", sep=",",head=TRUE, na.string="na")
#str(seed1)
###TRANFORMATION
skewness(seed1$weight)
seed_log = log1p(seed1$weight)
seed1 <- cbind(seed1,seed_log)

#DISTRIBUTION
hist(seed1$weight,breaks = 200)

#SUMMARIZE 
exp_sum <- summaryBy(pop+weight~pop+ssp+garden+year, data = seed1,FUN=c(mean,length))
exp_sum2 <- summaryBy(weight~ssp+garden+year, data = seed1,FUN=c(mean,min,max))

#write.csv(exp_sum2,file = "seed_exp_summary.csv")
#=================================================================================================#

###Boxplot
rr <- ggplot(seed1, aes(ssp,weight))
rr+geom_boxplot()+ facet_grid(.~garden)+ylim(-1,60)+xlab("subspecies")+theme_bw()

#Barplot (Fig.1)
ss <- ggplot(exp_sum2, aes(x=ssp,y=weight.mean,fill=factor(year, levels=c("y13","y12" ))))
ss + geom_bar(stat = "summary",fun.y = "mean")+facet_grid(.~garden)+xlab("subspecies")+
  ylab("seed weight")+ theme_bw(base_size = 15)+guides(fill=guide_legend(title=NULL))


#seedweight_means <- summaryBy(weight~ssp+garden, data= seed1, FUN = c(mean,sd,length))
#===================================================================================================#
###VARIABLE ELIMINATION
#POP MEANS
popgrw <- summaryBy(weight+MAT+MWMT+MCMT+TD+MAP+MSP+AHM+SHM+DD0+DD5+DDL18+DDG18+NFFD+bFFP+eFFP+FFP+
                      PAS+EMT+EXT+MAR+Eref+CMD+RH+PPT_wt+PPT_sp+PPT_sm+PPT_at+DD0_wt+DD5_sm ~ pop+garden+year+ssp+ssp, FUN = c(mean), data=seed1)
setnames(popgrw, old = c('weight.mean','MAT.mean','MWMT.mean','MCMT.mean','TD.mean','MAP.mean','MSP.mean','AHM.mean','SHM.mean','DD0.mean','DD5.mean','DDL18.mean','DDG18.mean','NFFD.mean','bFFP.mean','eFFP.mean','FFP.mean','PAS.mean','EMT.mean','EXT.mean','MAR.mean','Eref.mean','CMD.mean','RH.mean','PPT_wt.mean','PPT_sp.mean','PPT_sm.mean','PPT_at.mean','DD0_wt.mean','DD5_sm.mean'),
         new = c('weight','MAT','MWMT','MCMT','TD','MAP','MSP','AHM','SHM','DD0','DD5','DDL18','DDG18','NFFD','bFFP','eFFP','FFP','PAS','EMT','EXT','MAR','Eref','CMD','RH','PPT_wt','PPT_sp','PPT_sm','PPT_at','DD0_wt','DD5_sm'))

pop_sum <- summaryBy(weight+MAT+MWMT+MCMT+TD+MAP+MSP+AHM+SHM+DD0+DD5+DDL18+DDG18+NFFD+bFFP+eFFP+FFP+
                          PAS+EMT+EXT+MAR+Eref+CMD+RH+PPT_wt+PPT_sp+PPT_sm+PPT_at+DD0_wt+DD5_sm ~ pop+ssp+ssp, FUN = c(mean), data=seed1)

#write.xlsx(x = popgrw_sum,file = "data_pop.xlsx")

data_pop=with(popgrw, data.frame(MAT,MWMT,MCMT,TD,MAP,MSP,AHM,SHM,DD0,DD5,DDL18,DDG18,NFFD,bFFP,eFFP,FFP,PAS,EMT,EXT,MAR,Eref,CMD,RH,PPT_wt,PPT_sp,PPT_sm,PPT_at,DD0_wt,DD5_sm,weight)) 
cor <- cor(data_pop)

### HIGH VAR = DD0_wt, PPT_wt, no interaction
#=================================================================================================#
###GLMM 
sseed1 <- glmmTMB(
  weight ~ ssp + PPT_sm + (1 | garden:year) + (1|pop:(garden:year)), 
  ziformula = ~1, 
  family = tweedie(), 
  data = seed1)

sseed2 <- glmmTMB(
  weight ~ ssp + DD0_wt + PPT_sm + (1 | garden:year) + (1|pop:(garden:year)), 
  ziformula = ~1, 
  family = tweedie(), 
  data = seed1)

sseed3 <- glmmTMB(
  weight ~ ssp + (1 | garden:year) + (1|pop:(garden:year)), 
  ziformula = ~1, 
  family = tweedie(), 
  data = seed1)

#LRT for two models
anova(sseed1,sseed2,sseed3)
#BEST MODEL = sseed2
summary(sseed2)

###EXAMINE RESIDUALS
residuals <- resid(sseed2) 
summary(residuals)
hist(residuals)
plot(residuals)
testUniformity(sseed2)

#TABLE 2b
#With Climate
tab_model(sseed2, p.val = "wald", show.icc = FALSE, collapse.ci = TRUE)
#Without climate
tab_model(sseed3, p.val = "wald", show.icc = FALSE, collapse.ci = TRUE)
#GARDEN AND YEAR EFFECTS (SJPLOT)
#FigS1
plot_model(sseed1,type="re")

##GLMM OUTPUT
seed1_output <- augment(sseed2)
#CHANGE TO DATAFRAME
seed1_output <- as.data.frame(seed1_output)

#CONFIDENCE INTERVALS
confint(sseed2,level = 0.95,method = "wald")

#extract random effects
re1 <- ranef(sseed2)
popXgard <- (re1$cond$'pop:(garden:year)')

###SUMMARIZE FOR GxE
#BIND RANDOM EFFECTS 
seed_summary <- summaryBy(.fitted + weight + DD0_wt + PPT_sm ~ pop + ssp + garden +year,
                          data = seed1_output, FUN = c(mean))

#BIND CLIMATE VARS and RANDOM EFFECTS
seed_summary <- cbind(seed_summary,popXgard)
setnames(seed_summary, old = c('(Intercept)'),
         new = c('raneff'))
#SUMMARIZE WO YEAR
seed_summary1 <- summaryBy(weight.mean + .fitted.mean + raneff + DD0_wt.mean +
                             PPT_sm.mean  ~ pop + ssp + garden,data = seed_summary, FUN = c(mean))

#Relationship between PPT_sm and seed yield at each garden
ww <- ggplot(seed_summary1, aes(PPT_sm.mean.mean,weight.mean.mean,color=ssp,shape=ssp))
ww+ geom_point()+stat_smooth(method=lm,se=FALSE,linetype=4,fill=NA)+theme_bw()+
  xlab("Summer precipitation (mm)")+ylab("Predicted seed yield")+facet_grid(.~garden)

yy <- ggplot(seed_summary1, aes(DD0_wt.mean.mean,weight.mean.mean,color=ssp,shape=ssp))
yy+ geom_point()+stat_smooth(method=lm,se=FALSE,linetype=4,fill=NA)+theme_bw()+
  xlab("Winter degree-days < 0C")+ylab("Predicted seed yield")+facet_grid(.~garden)


#SUMMARIZE TO GARDEN EFFECTS 
Ogrw_GxE <- dplyr::filter(seed_summary1,garden=="Orchard")
Mgrw_GxE <- dplyr::filter(seed_summary1,garden=="Majors")

#GARDEN-PLOIDY CORRELATION
with(subset(Ogrw_GxE, ssp==c("t")),cor.test(PPT_sm.mean.mean, weight.mean.mean))
with(subset(Ogrw_GxE, ssp==c("w")),cor.test(PPT_sm.mean.mean, weight.mean.mean))
with(subset(Mgrw_GxE, ssp==c("t")),cor.test(PPT_sm.mean.mean, weight.mean.mean))
with(subset(Mgrw_GxE, ssp==c("w")),cor.test(PPT_sm.mean.mean, weight.mean.mean))
with(subset(Ogrw_GxE, ssp==c("t")),cor.test(DD0_wt.mean.mean, weight.mean.mean))
with(subset(Ogrw_GxE, ssp==c("w")),cor.test(DD0_wt.mean.mean, weight.mean.mean))
with(subset(Mgrw_GxE, ssp==c("t")),cor.test(DD0_wt.mean.mean, weight.mean.mean))
with(subset(Mgrw_GxE, ssp==c("w")),cor.test(DD0_wt.mean.mean, weight.mean.mean))

#PPT_sm RELATIONSHIP, SUBSET TRID at ORCHARD, SINCE SEED YIELD NEAR 0
tgrw_GxE <- dplyr::filter(seed_summary1,ssp=="t")

#Fig5 Relationship between PPT_sm and tridentata seed yield at each garden
ww <- ggplot(tgrw_GxE, aes(PPT_sm.mean.mean,weight.mean.mean,shape=ssp))
ww+ geom_point()+stat_smooth(method=lm,se=FALSE,linetype=4,fill=NA,color="gray")+theme_bw(base_size = 15)+
  labs(x= "Summer Precipitation", y=(bquote('Weight' ~(grams^-y)))) +facet_grid(.~garden)#+geom_text(aes(label=pop))

yy <- ggplot(tgrw_GxE, aes(DD0_wt.mean.mean,weight.mean.mean,shape=ssp))
yy+ geom_point()+stat_smooth(method=lm,se=FALSE,linetype=4,fill=NA,color="gray")+theme_bw(base_size = 15)+
  xlab("DD0_wt")+ylab("Weight")+facet_grid(.~garden)


#COR.TEST for PPT_sm at Orch
cor.test(tgrw_GxE$weight.mean.mean,tgrw_GxE$PPT_sm.mean.mean)



