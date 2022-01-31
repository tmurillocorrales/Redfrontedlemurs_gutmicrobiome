# Linear mixed model to test the drivers of bacteria alpha diversity for samples having fecal glucorticoid measurements (fgc) ####
#These scripts were developed with the support of Dr. Roger Mundry

#Load packages ####
library(car) 
library(lme4) 
library(MuMIn)

#Read metadata with Faith's Phylogenetic Diversity Index
xdata = read.delim("metadata_alphadiv_socint_para_feed.txt", header = T, sep = "\t")

#Remove all samples without fecal glucocorticoid measurements
test.data = xdata[xdata$n11oxo_CM_wet_feces != "no",]
#From 799 to 641 samples
#Convert column with fecal glucocorticoid measurements to numeric
test.data$n11oxo_CM_wet_feces <- as.numeric(as.character(test.data$n11oxo_CM_wet_feces)) 
#Confirm absence of repeated values
unique(test.data$n11oxo_CM_wet_feces)
str(test.data)

###### Check distributions for all the covariates ######
#Create vector with covariates of interest for plotting together
test.vars= c("ffr.prop", "fle.prop","rain", "ffl.prop","age_months", "soc.int", "n11oxo_CM_wet_feces", "richness.para")

par(mfrow= c(3,3))  #Split the windows for plotting

#Loop for creating scatter plots.
for(i in 1: length(test.vars)) {
  hist(test.data[,test.vars[i]], main =test.vars[i])}

## Some of them are a bit skewed proceed to log transform to attempt more symmetrical distributions.

par(mfrow= c(1,1)) #Return window to regular setup

##Age in months
#Check for zeros
sum(xdata$age_months==0)
#No zeros
#Compare distributions of original and log transform data
hist(xdata$age_months)
hist(log(xdata$age_months))
#Log transform and create new variable in dataframe
xdata$log.age_months = log(xdata$age_months)

##fgc
#Check for zeros
sum(test.data$n11oxo_CM_wet_feces==0)
#No zeros
#Compare distributions of original and log transform data
hist(test.data$n11oxo_CM_wet_feces)
hist(log(test.data$n11oxo_CM_wet_feces))
#Log transform and create new variable in dataframe
test.data$log.fgc = log(test.data$n11oxo_CM_wet_feces)

#The distribution of feeding proportions, social interactions, parasite richness and precipitation is similar after log transform, so leave the 
#covariates as they are.

#Check number of data points per categorical variable -> look for outliers & number of observations.
#Data points per individual
plot(table(test.data$individual))
#Data points per sex
plot(table(test.data$sex))
#Data points per age_category
plot(table(test.data$age_category))
#Few data points for juveniles & infants
plot(table(test.data$group))

#######Correlation#####

#Check for correlations between covariates through scatterplots
plot(test.data$fle.prop, test.data$ffr.prop)
plot(test.data$rain, test.data$fle.prop)
plot(test.data$rain, test.data$ffr.prop)
plot(test.data$fle.prop, test.data$ffl.prop)
plot(test.data$rain, test.data$ffl.prop)
plot(test.data$soc.int, test.data$ffr.prop)
plot(test.data$soc.int, test.data$fle.prop)
plot(test.data$soc.int, test.data$ffl.prop)
plot(test.data$soc.int, test.data$rain)
plot(test.data$richness.para, test.data$fle.prop)
plot(test.data$richness.para, test.data$ffl.prop)
plot(test.data$richness.para, test.data$ffr.prop)
plot(test.data$richness.para, test.data$rain)
plot(test.data$richness.para, test.data$log.fgc)

#Check for correlations between covariates using Pearson correlation
#Check for NAs in case there are no observations for certain individuals per month.
#Create vector with covariates to test.
test.vars = c("ffr.prop", "fle.prop","rain", "ffl.prop", "age_months", "soc.int", "n11oxo_CM_wet_feces", "richness.para")
#Calculate correlation
res=cor(test.data[,test.vars], method = "pearson", use = "complete.obs")
#Get the largest absolute correlation
max(abs(res[upper.tri(res)]))
#Round to three decimals.
round(res, 3)
#                    ffr.prop fle.prop   rain ffl.prop age_months soc.int n11oxo_CM_wet_feces richness.para
#ffr.prop               1.000   -0.133  0.029   -0.077     -0.005  -0.165              -0.058        -0.154
#fle.prop              -0.133    1.000 -0.013    0.021     -0.083  -0.056              -0.077         0.118
#rain                   0.029   -0.013  1.000    0.097      0.007  -0.064              -0.217        -0.054
#ffl.prop              -0.077    0.021  0.097    1.000      0.000  -0.076              -0.154         0.027
#age_months            -0.005   -0.083  0.007    0.000      1.000  -0.178              -0.024         0.116
#soc.int               -0.165   -0.056 -0.064   -0.076     -0.178   1.000               0.095        -0.150
#n11oxo_CM_wet_feces   -0.058   -0.077 -0.217*  -0.154     -0.024   0.095               1.000        -0.055
#richness.para         -0.154    0.118 -0.054    0.027      0.116  -0.150              -0.055         1.000

#No variables correlate as the highest value obtained is 0.2165882.

##### Random slopes ####
#Include random effect of individual and SampleID
#Test random slopes for feeding proportions, social interactions, parasite richness, age, fgc and precipitation

#Random slopes were determined by using a function developed by Dr. Roger Mundry.

source("diagnostic_fcns.r")

#Random slopes should be considered differently for factors & covariates.
#For a factor: include those with 2  or more observations per level of random effect.
#For a covariate: include those 3 or more unique observations of the fixed effect per level of the random effect. 
#IF the covariate does not meet the rule as a covariate, check if it meets it for a factor.
#If only some levels meet the requirement Roger recommends to include the random slope when at least half of the levels from the random effect do fulfill the rules.

#log.affr_rate_within_individual (covariate)`
#9  10  11  12 tot -> number of unique observations (here should be 3 or more)
#1   1   2   1   5 -> levels of the random effect

#log.affr_rate_within_individual (covariate as factor)`
#0   1 tot -> number of unique observations (here should be 2 or more)
#2   3   5 -> levels of the random effect

#Model for checking random slopes using function fe.re.tab
xx.fe.re = fe.re.tab(
  fe.model = "alpha.div ~ sex + group + log.age_months + soc.int + rain + fle.prop + ffr.prop + ffl.prop + log.fgc + richness.para",
  re = "(1|individual) + (1|SampleID)",
  data = test.data, other.vars = c("age_months", "n11oxo_CM_wet_feces")) 

xx.fe.re$summary

#According to the results the random slopes of feeding proportions, social interactions, parasite richness, age, and precipitationin individual should be included.
#Sex and group do not have random slopes.

#z transform the covariates to ease model converge.
t.data$z.fle.prop=as.vector(scale(t.data$fle.prop))
t.data$z.rain = as.vector(scale(t.data$rain))
t.data$z.ffr.prop = as.vector(scale(t.data$ffr.prop))
t.data$z.ffl.prop = as.vector(scale(t.data$ffl.prop))
t.data$z.log.age_months = as.vector(scale(t.data$log.age_months))
t.data$z.soc.int = as.vector(scale(t.data$soc.int))
t.data$z.log.fgc = as.vector(scale(t.data$log.fgc))
t.data$z.richness.para = as.vector(scale(t.data$richness.para))
str(t.data)

#Check distribution of the response variable alpha diversity (alpha.div)
hist(xdata$alpha.div)
#Log transform 
hist(log(xdata$alpha.div))
#After log transformation the distribution is not symmetrical.

#Transformation of PD ->Box-Cox power transformation
xx=lm(alpha.div ~ sex + group + z.log.age_months + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.soc.int + z.log.fgc + z.richness.para,
      data=t.data)

#Box.Cox transformation#####
xx=powerTransform(object=xx)
xx
#Applied the transformation such that we get the transformed response.
t.data$tr.alpha.div=(t.data$alpha.div^xx$lambda-1)/xx$lambda
transform.para=xx
#Rescale tr.PD to the same max as PD.
t.data$tr.alpha.div=t.data$tr.alpha.div/max(t.data$tr.alpha.div)
t.data$tr.alpha.div=t.data$tr.alpha.div*max(t.data$alpha.div)
plot(t.data$alpha.div, t.data$tr.alpha.div)

#To ease model converge choose the optimizer and increase optCtrl.
contr=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))

#Model with correlation ####
full.wc = lmer(tr.alpha.div ~ sex + group + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc + z.richness.para +
                 (1 + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc + z.richness.para|individual), REML=F,
               data=t.data, control = contr)

diagnostics.plot(full.wc)

#Residuals against fitted values are not perfect but still acceptable.

#Revise results from the model
summary(full.wc)$varcor

#Model without correlations #####

full.woc = lmer(tr.alpha.div ~ sex + group + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc + z.richness.para +
                  (1 + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc + z.richness.para||individual), REML=F,
                data=t.data, control = contr)

#Compare log likelihoods of both models to determine if correlations should be kept.
logLik(full.wc)#-2041.781 (df=59)
logLik(full.woc) #-2050.398 (df=23)
#Keep correlations as there is a big difference between both models.

#Collinearity ####
#Check for collinearity between predictors, this happens when values are higher than 5.
vif(full.wc)
vif(full.wc)[,3]^2
#No collinearity as the highest value was 1.404524 for parasite richness.

#Check distribution of random effects -> BLUPs (Best linear unbiased predictors).
ranef.diagn.plot(full.wc)
#They are acceptable.

# Model stability ####
#Model stability was determined by using a function developed by Dr. Roger Mundry.
source("glmm_stability.r")

#Function for checking model stability
full.stab = glmm.model.stab(model.res = full.wc, contr = NULL, para = F, data = NULL)

#Evaluate warnings 
table(full.stab$detailed$warnings)
#They are acceptable

#Evaluate stability
round(full.stab$summary [,-1], 3)

#Prepare graphical representation for better understanding.
is.re = grepl(x=row.names(full.stab$summary), pattern = "@")
is.re

#Stability of fixed effects:
m.stab.plot(full.stab$summary[!is.re, -1])
#Acceptable

#Stability of random effects:
m.stab.plot(full.stab$summary[is.re, -1])
##They are acceptable with the exception of the individual.

# Inference ##############
summary(full.wc)

#Full-null model comparison #####

###Null model ####
#Prepare a null model excluding only the predictors of interest group, social interactions, age, sex, fgc and parasite richness from 
#the fixed effects but including their random slopes.

#Change optimizer to ease model converge
contr=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))

null = lmer(tr.alpha.div ~ z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop +
              (1 + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc + z.richness.para|individual), REML=F,
            data=t.data, control = contr)


#Comparison between full-null model.
as.data.frame(anova(null, full.wc, test="Chisq"))
#pvalue =  0.03890784
#The full and the null models are significantly different the predictors of interest impact alpha diversity.
#Results shown in supplementary table S9

# Test the impact of single predictors ####
#lmertest ####
library(lmerTest) #-> how to test single predictors in Gaussian models.

#Set optmizer to ease model convergence
contr.reml=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)) #optimizers nlminbwrap bobyqa

#Calculate a new model in lmerTest for testing the effect of single predictors.
full.wc.reml = lmer(tr.alpha.div ~ sex + group + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc + z.richness.para +
                      (1 + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc + z.richness.para|individual), REML=F,
                    data=t.data, control = contr.reml)

#Estimates full model in lmer
round(summary(full.wc)$coefficients,3)
#Estimates full model in lmerTest
round(summary(full.wc.reml)$coefficients,3)
#Even though it did not converge the estimates are very similar to the original model calculated with lmer, so the results from lmertest can be used.

#                   Estimate Std. Error      df t value Pr(>|t|)
#(Intercept)        26.786      1.506  31.021  17.782    0.000
#sexmale             0.801      0.997 197.154   0.804    0.423
#groupB              0.074      1.728  31.751   0.043    0.966
#groupF             -5.210      1.907  53.987  -2.733    0.008*
#groupJ             -1.777      1.711  34.705  -1.038    0.306
#z.log.age_months   -0.001      0.510 204.977  -0.002    0.999
#z.soc.int           0.613      0.618   8.713   0.992    0.348
#z.rain              0.751      0.720  19.765   1.043    0.310
#z.fle.prop          2.205      0.556  28.125   3.969    0.000*
#z.ffr.prop         -0.440      0.580  27.561  -0.757    0.455
#z.ffl.prop         -0.700      0.560   6.979  -1.250    0.252
#z.log.fgc           1.216      0.532  45.259   2.285    0.027*
#z.richness.para    -0.871      0.653  28.834  -1.334    0.193

#Export data
estimates = as.data.frame(round(summary(full.wc.reml)$coefficients,3))
#Export table for results shown in supplementary table S9
write.table(estimates, file = "estimates_alpha_para_fgc.csv", sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)

detach(name = "package:lmerTest") #important when lmer will be used again.

# Confidence intervals #####
#Confidence intervals were calculated by using a function developed by Dr. Roger Mundry.
source("boot_glmm.r")

#ci.estimates confidence limits from estimates
#ci.predicted confidence limits of the fitted model over a range of values of the predictors
#all.warns warnings
#all.boots results of individual bootstraps
#If warning "closing unused connection" comes up it can be ignored.

#Package needed for the bootstrapping function
library(kyotil)

#Calculation of confidence intervals required high computational capacities.
#Calculate model
boot.full = boot.lmer(m = full.wc, discard.warnings = F, nboots = 1000, para = T, n.cores = 5,
          resol = 1000, level = 0.95, use = "group")

#Check results
round(boot.full$ci.estimates, 3)
#Plot results
m.stab.plot(boot.full$ci.estimates)
#Check results
boot.full$ci.fitted
#Confidence intervals are acceptable

# Plots #####
#Plotting was done by using two functions developed by Dr. Roger Mundry.
source("factor_int_plot_fcn.r")
source("no_int_plot.r")

t.data$group.B= t.data$group.B-mean(t.data$group.B)
t.data$group.F= t.data$group.F-mean(t.data$group.F)
t.data$group.J= t.data$group.J-mean(t.data$group.J)
t.data$sex.male = t.data$sex.male-mean(t.data$sex.male)
names(t.data)

#Calculate model for plot
plot.full.wc = lmer(tr.alpha.div ~ sex.male + group.B + group.F + group.J + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc+ z.richness.para +
                      (1 + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.log.fgc + z.richness.para|individual), REML=F,
                    data=t.data, control = contr)

#Extract confidence intervals for feeding on leaves and fecal glucocorticoids
plot.boot.full = boot.lmer(m = plot.full.wc, discard.warnings = F, nboots = 1000, para = T, n.cores = 5,
                           resol = 1000, level = 0.95, use = list("z.fle.prop", "z.log.fgc"))  


range(t.data$alpha.div)
yaxt.labels=seq(30,50, by=10)
yaxt.at=(yaxt.labels^transform.para$lambda-1)/transform.para$lambda
yaxt.at=yaxt.at/max((t.data$alpha.div^transform.para$lambda-1)/transform.para$lambda)
yaxt.at=yaxt.at*max(t.data$alpha.div)


##Plot fecal glucocorticoids -> figure 2c
#Set confidence intervals
ci=plot.boot.full$ci.fitted$z.log.fgc
#Set axis
range(t.data$n11oxo_CM_wet_feces)
x.labels=c(10,50, 100, 500, 1000)
x.at=(log(x.labels)-mean(log(t.data$n11oxo_CM_wet_feces)))/sd(log(t.data$n11oxo_CM_wet_feces))
#Plot
no.int.plot (
  plot.data = t.data, covariate = "z.log.fgc", coefs=fixef(plot.full.wc), response = t.data$tr.alpha.div, link=c("identity"),
  grid.resol=NA, 
  xlab="Fecal glucocorticoids", ylab="Alpha diversity", x.at=x.at, x.labels=x.labels, y.at=yaxt.at, y.labels=yaxt.labels, xlim=NULL, ylim=NULL,
  cex.lab=1, cex.axis=1, size.fac=1.5, col=grey(level=0.5, alpha = 0.5), pch=19, conf.int=ci, ci.type=c("lines", "area"), lty=2, lwd=1, pt.lwd=1,
  my.par=list(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15, las=1), quiet=T, reset.par=F, mean.fun=mean,
  weights=NULL)