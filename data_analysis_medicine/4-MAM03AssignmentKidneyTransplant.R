library(survival)  # core survival analysis function
library(survminer)
library(splines)
library(dplyr)
library(tidyr)
library(table1) 
library(car)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(rms)

setwd("/home/dilaratank/Desktop/Vakken/MSc/Year 1/Advanced data analysis in medicine/week 4")
load("kidney transplant.RData")

# COX REGRESSION: combine events into one
d$stat_pat = as.numeric(d$stat_pat)-1
d$stat_gra = as.numeric(d$stat_gra)-1

d$dfs = d$time_to_graft_failure
d$dfs[d$stat_gra==0] = d$ time_to_death[d$stat_gra==0]
d$dfs_stat = d$stat_gra
d$dfs_stat[d$stat_pat==1] = 1

# so that you can use: coxph(Surv(dfs, dfs_stat) ~ agedon+sex_pat+screat1+ .... , data=d)

# COX REGRESSION: linearity check
covs <- data.frame(d[, c("bmi", "age_at_tx", "duur_dia", "agedon", "gfr3", "screat3", "map3")])

library(splines)
par(mfrow=c(4,2))
psplinespvals=c()
for (i in 1:ncol(covs)) {
  x=qq=NA
  x=as.numeric(unlist(covs[i]))         # select the covariate of interest
  if (is.numeric(x) ==TRUE) {
    aic=c()
    for (j in 1:12) {aic[j] = AIC(coxph(Surv(dfs, dfs_stat) ~ ns(x, df=j), data=d))}
    plot(aic, type="b", main=names(covs)[i], xlab="natural spline penalty (df)")
    qq=summary(coxph(Surv(dfs, dfs_stat) ~ pspline(x), data=d)) 
    psplinespvals[i] = qq$coefficients[2,6]
    text(x=1, y=min(aic,na.rm=T)+(max(aic,na.rm=T)-min(aic,na.rm=T))/2, 
         paste("pspline nonlin p = ", round(psplinespvals[i],4)), adj=0, col=4, cex=0.9)
  }
}

covs <- data.frame(d[, c("age_at_tx", "duur_dia", "gfr3")])

library(splines)
par(mfrow=c(3,1))
psplinespvals=c()
for (i in 1:ncol(covs)) {
  x=qq=NA
  x=as.numeric(unlist(covs[i]))         # select the covariate of interest
  if (is.numeric(x) ==TRUE) {
    aic=c()
    for (j in 1:12) {aic[j] = AIC(coxph(Surv(dfs, dfs_stat) ~ ns(x, df=j), data=d))}
    plot(aic, type="b", main=names(covs)[i], xlab="natural spline penalty (df)")
    qq=summary(coxph(Surv(dfs, dfs_stat) ~ pspline(x), data=d)) 
    print(qq)
    psplinespvals[i] = qq$coefficients[2,6]
    text(x=1, y=min(aic,na.rm=T)+(max(aic,na.rm=T)-min(aic,na.rm=T))/2, 
         paste("pspline nonlin p = ", round(psplinespvals[i],4)), adj=0, col=4, cex=0.9)
  }
}

# COX REGRESSION: proportionality check
label(d$sex_pat) <- "Sex of patient"
label(d$bmi) <- "BMI"
label(d$agedon) <- "Age of donor"
label(d$sexdon) <- "Sex of donor"
label(d$type_dia) <- "Type of dialysis"
label(d$duur_dia) <- "Duration of dialysis"
label(d$retrans) <- "Re-transplantation"
label(d$gfr1) <- "Glomerus filtration rate at 3 months"
label(d$gfr2) <- "Glomerus filtration rate at 6 months"
label(d$gfr3) <- "Glomerus filtration rate at 1 year"
label(d$gfr4) <- "Glomerus filtration rate at 2 years"
label(d$gfr5) <- "Glomerus filtration rate at 5 years"
label(d$gfr6) <- "Glomerus filtration rate at 10 years"
label(d$retrans) <- "Re-transplantation"
label(d$screat1) <- "Serum creatinine at 3 months"
label(d$screat2) <- "Serum creatinine at 6 months"
label(d$screat3) <- "Serum creatinine at 1 year"
label(d$screat4) <- "Serum creatinine at 2 years"
label(d$screat5) <- "Serum creatinine at 5 years"
label(d$screat6) <- "Serum creatinine at 10 years"
label(d$map1) <- "Mid arterial pressure at 3 months"
label(d$map2) <- "Mid arterial pressure at 6 months"
label(d$map3) <- "Mid arterial pressure at 1 year"
label(d$map4) <- "Mid arterial pressure at 2 years"
label(d$map5) <- "Mid arterial pressure at 5 years"
label(d$map6) <- "Mid arterial pressure at 10 years"
label(d$age_at_tx) <- "Age at transplantation"
label(d$time_to_graft_failure) <- "Years from transplantation till date of graft failure"
label(d$time_to_death) <- "Years from transplantation till date of death"

cox_multi <- (coxph(Surv(dfs,dfs_stat) ~ as.factor(sex_pat)+ns(age_at_tx,df=2)+as.factor(type_dia)+ns(duur_dia, df=2)+bmi+as.factor(sexdon)+agedon+as.factor(retrans)+ns(gfr3, df=9)+screat3+map3, data=d))
cox.zph(cox_multi)
cox_multi

labels = c("Sex of patient", 
           "Age at transplantation: 1", "Age at transplantation: 2", 
           "Type of dialysis: geen dialyse", "Type of dialysis: Hemodialyse centrum", "Type of dialysis: Hemodialyse thuis", 
           "Duration of dialysis: 1", "Duration of dialysis: 2", 
           "BMI", 
           "Sex of donor", 
           "Age of donor", 
           "Re-transplantation", 
           "Glomerus filtration rate at 1 year: 1", "Glomerus filtration rate at 1 year: 2", "Glomerus filtration rate at 1 year: 3", "Glomerus filtration rate at 1 year: 4", "Glomerus filtration rate at 1 year: 5", "Glomerus filtration rate at 1 year: 6", "Glomerus filtration rate at 1 year: 7", "Glomerus filtration rate at 1 year: 8", "Glomerus filtration rate at 1 year: 9", 
           "Serum creatinine at 1 year", 
           "Mid arterial pressure at 1 year")
tab_model(cox_multi, string.est = "Hazard Ratio", pred.labels=labels)
summary(cox_multi)
Cstat(cox_multi)

# LME: remove outliers gfr
dlong1=subset(dlong,is.na(gfr)==FALSE)      # verwijder records met missing gfr-metingen
dlong1$gfr[dlong1$gfr > 200]=NA             # verwijder 2 records met onwaarschijnlijke (onmogelijk) hoge gfr metingen
dlong1=subset(dlong,is.na(gfr)==FALSE)

# LME: gfr
lme_gfr=lme(gfr~ns(years,df=2),random=~1+ns(years,df=2)|ID,data=dlong1,method="REML")
summary(lme_gfr)

# LME: remove NAN map
dlong1=subset(dlong,is.na(map)==FALSE)      # verwijder records met missing map-metingen

# LME: map
lme_map=lme(map~ns(years,df=2),random=~1+ns(years,df=2)|ID,data=dlong1,method="REML",control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
summary(lme_map)

# LME: assumption check map
par(mfrow=c(3,2))
plot((dlong1$years), resid(lme_map),main="Constant variance of residuals",xlab="years follow-up",ylab="Residuals")
abline(h=0, lty=2)
lines(unique(round(dlong1$years))[order(unique(round(dlong1$years)))],
      tapply(resid(lme_map), round(dlong1$years), mean)[order(unique(round(dlong1$years)))],
      col=2, lwd=3)
qqnorm(resid(lme_map),main="Constant variance of residuals")
qqline(resid(lme_map))

hist(lme_map$coefficients$random$ID[,1],xlab="Patient specific intercept",main="Normality of random intercept")
qqnorm(lme_map$coefficients$random$ID[,1],main="Normality of random intercept")
qqline(lme_map$coefficients$random$ID[,1])
hist(lme_map$coefficients$random$ID[,2],xlab="Patient specific slope",main="Normality of random slope")
qqnorm(lme_map$coefficients$random$ID[,2],main="Normality of random slope")
qqline(lme_map$coefficients$random$ID[,2])


# LME: remove NAN screat
dlong1=subset(dlong,is.na(creat)==FALSE)      # verwijder records met missing map-metingen

# LME: screat
lme_creat=lme(ns(creat, df=2)~ns(years,df=2),random=~1+ns(years,df=2)|ID,data=dlong1,method="REML",control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
summary(lme_creat)

par(mfrow=c(3,2))
plot((dlong1$years), resid(lme_creat),main="Constant variance of residuals",xlab="years follow-up",ylab="Residuals")
abline(h=0, lty=2)
lines(unique(round(dlong1$years))[order(unique(round(dlong1$years)))],
      tapply(resid(lme_creat), round(dlong1$years), mean)[order(unique(round(dlong1$years)))],
      col=2, lwd=3)
qqnorm(resid(lme_creat),main="Constant variance of residuals")
qqline(resid(lme_creat))

hist(lme_creat$coefficients$random$ID[,1],xlab="Patient specific intercept",main="Normality of random intercept")
qqnorm(lme_creat$coefficients$random$ID[,1],main="Normality of random intercept")
qqline(lme_creat$coefficients$random$ID[,1])
hist(lme_creat$coefficients$random$ID[,2],xlab="Patient specific slope",main="Normality of random slope")
qqnorm(lme_creat$coefficients$random$ID[,2],main="Normality of random slope")
qqline(lme_creat$coefficients$random$ID[,2])


library(JMbayes2)

dlong1=subset(dlong,is.na(gfr)==FALSE & is.na(creat)==FALSE & is.na(map)==FALSE)  # verwijder alle records met missing

dlong1=subset(dlong1,(dlong1$years>=dlong1$time_to_death)==FALSE)

dlong1$gfr[dlong1$gfr > 200]=NA    # verwijder die 2 rare gfr metingen

dlong1=subset(dlong1,is.na(gfr)==FALSE)

dim(dlong1)   # check hoeveel records er in dlong1 zitten

lmemodel_gfr = lme(gfr ~ ns(years,2)+sex_pat+bmi+sexdon+agedon+type_dia+duur_dia+retrans+age_at_tx, random = ~1+ns(years,2)|ID, data=dlong1)  

summary(lmemodel_gfr)    # check hoeveel pats en hoeveel records in de analyse zitten (dit moet hetzelfde aantal als bij de 2 andere lme's zijn)

lmemodel_creat = lme(creat ~ years+sex_pat+bmi+sexdon+agedon+type_dia+duur_dia+retrans+age_at_tx, random = ~1+years|ID, data=dlong1)

summary(lmemodel_creat)

lmemodel_map = lme(map ~ years+sex_pat+bmi+sexdon+agedon+type_dia+duur_dia+retrans+age_at_tx, random = ~1+years|ID, data=dlong1)

summary(lmemodel_map)

jmresult = jm(cox_multi, list(lmemodel_gfr, lmemodel_creat, lmemodel_map), time_var="years")

summary(jmresult)


#JOINT MODEL
lmemodel_gfr = lme(gfr ~ ns(years,2), random = ~1+ns(years,2)|ID, data=dlong1)  

summary(lmemodel_gfr)    # check hoeveel pats en hoeveel records in de analyse zitten (dit moet hetzelfde aantal als bij de 2 andere lme's zijn)

lmemodel_creat = lme(creat ~ years, random = ~1+years|ID, data=dlong1)

summary(lmemodel_creat)

lmemodel_map = lme(map ~ years, random = ~1+years|ID, data=dlong1)

summary(lmemodel_map)

jmresult = jm(cox_multi, list(lmemodel_gfr, lmemodel_creat, lmemodel_map), time_var="years")

summary(jmresult)




