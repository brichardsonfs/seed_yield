require(lmerTest)
require(ggplot2)
require(doBy)
require(MuMIn)
require(lattice)
require(broom)
require(data.table)
require(sjPlot)
require(e1071)
require(insight)
require(dplyr)
require(lattice)
require(Hmisc)
#==============================================================================================#
###READ IN DATA
seed1 <- read.csv(file="seedmaster_climNA.csv", sep=",",head=TRUE, na.string="na")
#str(seed1)
###TRANFORMATION
skewness(seed1$weight)
seed_log = log1p(seed1$weight)
seed1 <- cbind(seed1,seed_log)


#SUM ACROSS YEARS (2011 to 2012)
#seed1yr <- summaryBy(weight + seed_log ~ ~ pop + ssp + ploidy + family + garden+MAT+MWMT+MCMT+
#                       TD+MAP+MSP+AHM+SHM+DD0+DD5+DDL18+DDG18+NFFD+bFFP+eFFP+FFP+PAS+EMT+EXT+
#                       MAR+Eref+CMD+RH+PPT_wt+PPT_sp+PPT_sm+PPT_at+DD0_wt+DD5_sm,data = seed1,FUN = c(sum))
#setnames(seed1yr, old = c('weight.sum','seed_log.sum'), new = c('weight','seed_log'))
#SUMMARIZE 
exp_sum <- summaryBy(pop+weight~pop+ssp+garden+year, data = seed1,FUN=c(mean,length))
exp_sum2 <- summaryBy(weight~ssp+garden+year, data = seed1,FUN=c(mean,min,max))

write.csv(exp_sum2,file = "seed_exp_summary.csv")
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
popgrw <- summaryBy(seed_log+MAT+MWMT+MCMT+TD+MAP+MSP+AHM+SHM+DD0+DD5+DDL18+DDG18+NFFD+bFFP+eFFP+FFP+
                      PAS+EMT+EXT+MAR+Eref+CMD+RH+PPT_wt+PPT_sp+PPT_sm+PPT_at+DD0_wt+DD5_sm ~ pop+garden+year+ssp+ssp, FUN = c(mean), data=seed1)
setnames(popgrw, old = c('seed_log.mean','MAT.mean','MWMT.mean','MCMT.mean','TD.mean','MAP.mean','MSP.mean','AHM.mean','SHM.mean','DD0.mean','DD5.mean','DDL18.mean','DDG18.mean','NFFD.mean','bFFP.mean','eFFP.mean','FFP.mean','PAS.mean','EMT.mean','EXT.mean','MAR.mean','Eref.mean','CMD.mean','RH.mean','PPT_wt.mean','PPT_sp.mean','PPT_sm.mean','PPT_at.mean','DD0_wt.mean','DD5_sm.mean'),
         new = c('seed_log','MAT','MWMT','MCMT','TD','MAP','MSP','AHM','SHM','DD0','DD5','DDL18','DDG18','NFFD','bFFP','eFFP','FFP','PAS','EMT','EXT','MAR','Eref','CMD','RH','PPT_wt','PPT_sp','PPT_sm','PPT_at','DD0_wt','DD5_sm'))

pop_sum <- summaryBy(seed_log+MAT+MWMT+MCMT+TD+MAP+MSP+AHM+SHM+DD0+DD5+DDL18+DDG18+NFFD+bFFP+eFFP+FFP+
                          PAS+EMT+EXT+MAR+Eref+CMD+RH+PPT_wt+PPT_sp+PPT_sm+PPT_at+DD0_wt+DD5_sm ~ pop+ssp+ssp, FUN = c(mean), data=seed1)

#write.xlsx(x = popgrw_sum,file = "data_pop.xlsx")

data_pop=with(popgrw, data.frame(MAT,MWMT,MCMT,TD,MAP,MSP,AHM,SHM,DD0,DD5,DDL18,DDG18,NFFD,bFFP,eFFP,FFP,PAS,EMT,EXT,MAR,Eref,CMD,RH,PPT_wt,PPT_sp,PPT_sm,PPT_at,DD0_wt,DD5_sm,seed_log)) 
cor <- cor(data_pop)

### VAR w/ > +/- 0.3 = DD0,EMT, MCMT, PPT_wt
## MCMT and PPT_wt highest non-colinear vars

#data_a=with(seed1yr, data.frame(MAT,MWMT,MCMT,TD,MAP,MSP,AHM,SHM,DD0,DD5,DDL18,DDG18,NFFD,bFFP,eFFP,FFP,PAS,EMT,EXT,MAR,Eref,CMD,RH,PPT_wt,PPT_sp,PPT_sm,PPT_at,DD0_wt,DD5_sm,seed_log))

###OMIT VARIABLES < -/+0.25; DD0_wt, MCMT, EMT, DD0, PPT_wt; due to collinearity >0.65 EMT AND PPT_wt AND PPT_sm REMAIN
full1 <- lmer (seed_log ~ + PPT_sm + EMT + PPT_wt + (1|year:garden) + (1|pop:(garden:year)),REML=TRUE , data=seed1)
s_full <- step(full1) 
s_full


##GENETIC MODEL
aseed1 <- lmer (seed_log ~  ssp + MCMT + (1|garden:year) + (1|pop:(garden:year)), REML=TRUE , data=seed1)
aseed2 <- lmer (seed_log ~  ssp + PPT_sm + (1|garden:year) + (1|pop:(garden:year)), REML=TRUE , data=seed1)
aseed3 <- lmer (seed_log ~  ssp + MCMT + PPT_sm + (1|garden:year) + (1|pop:(garden:year)), REML=TRUE , data=seed1)
aseed_wo_clim <- lmer (seed_log ~  ssp + (1|garden:year) + (1|pop:(garden:year)), REML=TRUE , data=seed1)

###TEST MODELS
anova(aseed1,aseed2,aseed3)

##CHOSE ASEED3 BASED ON AIC
##GENETIC MODEL
summary(aseed3)
ranova(aseed3)
r.squaredGLMM(aseed3)
# VAR EXPLAINED WITHOUT CLIMATE
r.squaredGLMM(aseed_wo_clim)
aseed3
#Change to lmerMod for tables
class(aseed3) <- "lmerMod"
#TABLE 2b
tab_model(aseed3, p.val = "kr", show.icc = FALSE)

##LMER OUTPUT
seed1_output <- augment(aseed1)

###EXAMINE RESIDUALS
residuals <- resid(aseed1) 
summary(residuals)
hist(residuals)
plot(residuals)

###EXTRACT FIXED EFFECTS
#y.hat4 <- model.matrix(aseed3 , ssp = "fixed") %*% fixef(aseed3)
#y.hat4ssp <- model.matrix(aseed0 , ssp = "fixed") %*% fixef(aseed0)

###EXTRACT RANDOM EFFECTS
re_pop <- ranef(aseed1,condVar=TRUE, whichel = "pop:(garden:year)")
re_env <- ranef(aseed1,condVar=TRUE, whichel = "garden:year")
#Fig.4
dotplot(re_env)
#dotplot(re_pop)

#EXAM POP INTERACTION
re_pop1 <- unlist(re_pop)
re_pop2 <- as.vector(re_pop1)

###SUMMARIZE FOR GxE
#BIND RANDOM EFFECTS 
seed_summary <- summaryBy(.fitted + .fixed + seed_log ~ pop + ssp + garden +year,
                          data = seed1_output, FUN = c(mean))

#BIND CLIMATE VARS and RANDOM EFFECTS
seed_summary <- cbind(seed_summary,popgrw$MCMT,popgrw$PPT_sm,re_pop2)
setnames(seed_summary, old = c('popgrw$MCMT','popgrw$PPT_sm'),
         new = c('MCMT','PPT_sm'))
#SUMMARIZE WO YEAR
seed_summary1 <- summaryBy(seed_log.mean + re_pop2+ MCMT +
                             PPT_sm ~ pop + ssp + garden,data = seed_summary, FUN = c(mean))

#Relationship between MCMT and seed yield at each garden
#ww <- ggplot(seed_summary1, aes(MCMT.mean,seed_log.mean.mean,color=ssp,shape=ssp))
#ww+ geom_point()+stat_smooth(method=lm,se=FALSE,linetype=4,fill=NA)+theme_bw()+
#  xlab("MCMT")+ylab("log(Seed Yield)")+facet_grid(.~garden)

#SUMMARIZE TO GARDEN EFFECTS 
Ogrw_GxE <- dplyr::filter(seed_summary1,garden=="Orchard")
Mgrw_GxE <- dplyr::filter(seed_summary1,garden=="Majors")

#GARDEN-SSP CORRELATION
with(subset(Ogrw_GxE, ssp==c("t")),cor.test(PPT_sm.mean, seed_log.mean.mean))
with(subset(Ogrw_GxE, ssp==c("w")),cor.test(PPT_sm.mean, seed_log.mean.mean))
with(subset(Mgrw_GxE, ssp==c("t")),cor.test(PPT_sm.mean, seed_log.mean.mean))
with(subset(Mgrw_GxE, ssp==c("w")),cor.test(PPT_sm.mean, seed_log.mean.mean))
with(subset(Ogrw_GxE, ssp==c("t")),cor.test(MCMT.mean, seed_log.mean.mean))
with(subset(Ogrw_GxE, ssp==c("w")),cor.test(MCMT.mean, seed_log.mean.mean))
with(subset(Mgrw_GxE, ssp==c("t")),cor.test(MCMT.mean, seed_log.mean.mean))
with(subset(Mgrw_GxE, ssp==c("w")),cor.test(MCMT.mean, seed_log.mean.mean))

#PPT_sm RELATIONSHIP, SUBSET TRID at ORCHARD, SINCE SEED YIELD NEAR 0
tgrw_GxE <- dplyr::filter(seed_summary1,ssp=="t")

#Fig5 Relationship between PPT_sm and tridentata seed yield at each garden
ww <- ggplot(tgrw_GxE, aes(PPT_sm.mean,seed_log.mean.mean))
ww+ geom_point()+stat_smooth(method=lm,se=FALSE,linetype=4,fill=NA,color="gray")+theme_bw(base_size = 15)+
  xlab("Summer Precipitation")+ylab("log(Seed Yield)")+facet_grid(.~garden)#+geom_text(aes(label=pop))

yy <- ggplot(Mgrw_GxE, aes(MCMT.mean,seed_log.mean.mean))
yy+ geom_point()+stat_smooth(method=lm,se=FALSE,linetype=4,fill=NA,color="gray")+theme_bw(base_size = 15)+
  xlab("Mean Coldest Month Temperature")+ylab("log(Seed Yield)")+facet_grid(.~ssp)


#COR.TEST for PPT_sm at Orch
cor.test(seed_summaryOrch$seed_log.mean.mean,seed_summaryOrch$PPT_sm.mean)
#COR.TEST for TD at Majors Not Significant
cor.test(seed_summaryMaj$seed_log.mean.mean,seed_summaryMaj$TD.mean)


##MODEL FIT
nn <- ggplot(seed4_summary, aes(weight.mean,expm1(.fixed.mean)))
nn + theme_bw()+stat_smooth(method=lm,se=FALSE,linetype=4, color="gray")+ 
  geom_point(aes(shape=ssp))+labs(x="Observed", y="Predicted",check_overlap=TRUE)+
  scale_shape_manual(values=c(20,23))+facet_grid(.~garden)


fit <- with(seed1yr2, data.frame(pop=pop, garden=garden, ssp=ssp, ssp=ssp, family=family, EMT, PPT_sm, TD, observed=weight,fitted=fitted(aseed3)))
fit <- cbind(fit,y.hat4) 

fit_pop <- summaryBy(observed + fitted + y.hat4 + EMT + PPT_sm + TD ~ pop + ssp + garden, data= fit, FUN = c(mean))
fit_pop <- cbind(fit_pop, dd$condval,dd$condsd)
write.csv(x = fit_pop, file = "seed_lmer_data.csv")
###ADDED CLIMATE DISTANCE FROM SEED SOURCE TO GARDEN
fit_pop1 <- read.csv(file="seed_lmer_data1.csv",sep=",",head=TRUE)

####PLOT GxE
r <- ggplot(fit_pop,aes(TD.mean,dd$condval))+theme_bw()+ scale_shape(solid=FALSE)
r + facet_grid(garden~.) + geom_point(aes(shape=ssp)) + stat_smooth(method=lm,se=FALSE,linetype=4, color="gray")

###SUMMARIZE FOR SSP

fitssp <- with(seed1yr2, data.frame(pop=pop, garden=garden, ssp=ssp, ssp=ssp, family=family, observed=weight,fitted=fitted(aseed0)))
fitssp <- cbind(fitssp,y.hat4ssp) 

fit_pop_ssp <- summaryBy(observed + fitted + y.hat4ssp ~ pop + ssp + garden, data= fitssp, FUN = c(mean))
fit_pop_ssp <- cbind(fit_pop_ssp, ss$condval,ss$condsd)

###GARDEN EFFECTS AT SSP LEVEL
ee <- ggplot(fit_pop_ssp, aes(y=ss$condval,x=garden,group=pop))
ee + geom_point()+ geom_text(aes(label=pop))+ geom_line()+facet_grid(.~ssp) + theme_bw()

ee <- ggplot(seed1yr2, aes(EMT,seed_log, color=ssp))
ee + geom_point()



#geom_errorbarh((aes(xmin=condval -2*condsd,xmax=condval +2*condsd)))
###GARDEN-YEAR EFFECTS AT POP
re_pop <- summaryBy(condval+condsd~ssp+garden+pop,data = re,FUN = c(mean))

###GARDEN-YEAR EFFECTS AT SSP
re_ssp <- summaryBy(condval+condsd~ssp+garden+year,data = re,FUN = c(mean))

e <- ggplot(re_pop, aes(y=condval.mean,x=garden,group=pop))
e + geom_point()+ geom_text(aes(label=pop,hjust=1, vjust=0))+ geom_line()+facet_grid(.~ssp) + theme_bw()

ee <- ggplot(re_ssp, aes(x=year,y=condval.mean, color=ssp))
ee + geom_point()+ facet_grid(garden~.) + geom_pointrange((aes(ymin=condval.mean -2*condsd.mean,ymax=condval.mean +2*condsd.mean)))+theme_bw()




###CONFIDENCE INTERVALS
#pp <- profile(aseed3)
#ci_0.4 <- confint(pp,level=0.4)

#USING PREDICTIONINTERVAL from MERTOOLS

PI <- predictInterval(merMod = aseed3, newdata = seed1, level = 0.4, n.sims = 1000, stat = "median",
                      ssp="linear.prediction", include.resid.var = TRUE)
PIwpops <- cbind(seed1$pop,PI)
colnames(PIwpops)[1] <- "pop"
PIavg <- summaryBy(fit+ upr + lwr ~ pop, data = PIwpops, FUN = c(mean))

###CALCULATION of AVERAGE 2X CI
CI2X <- PIavg$upr.mean-PIavg$lwr.mean
mean(CI2X)
#==============================================================================================#

###PLOT OF CIs FOR 39 POPULATIONS
gg <-ggplot(aes(x=pop, y=fit.mean, ymin=lwr.mean, ymax=upr.mean), data=PIavg)
gg + geom_point() + geom_linerange() + labs(x="Populations", y="Prediction w/ 40% PI",check_overlap=TRUE) + theme_bw()

#==============================================================================================#




#==============================================================================================#

###FIXED EFF GRAPH FOR MANUSCRIPT

q<-ggplot(fit_pop_g, aes(y=observed.mean,x=fitted.mean))+theme_bw() 
q+stat_smooth(method=lm,se=FALSE,linessp=4,color="gray",size=1) + geom_point(aes(shape=garden),size=3) + xlab("Predicted") + ylab("Observed") + labs(shape = "Gardens") +scale_shape(solid=FALSE) 

p<-ggplot(fit_pop_g, aes(y=observed.mean,x=y.hat4.mean,shape=garden,fill=garden))+theme_bw() 
p+stat_smooth(method=lm,se=FALSE,linessp=4,size=1,color="gray") + geom_point(size=3) + xlab("Predicted") + ylab("Observed") + labs(color = "Gardens") + scale_shape_manual(values=c(21,22,24))
