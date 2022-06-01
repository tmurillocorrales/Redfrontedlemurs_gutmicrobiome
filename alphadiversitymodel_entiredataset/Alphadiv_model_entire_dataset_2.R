# Linear mixed model to test the drivers of bacteria alpha diversity ####
#These scripts were developed with the support of Dr. Roger Mundry

#Load packages ####
library(car) 
library(lme4) 
library(MuMIn)

#Read metadata with Faith's Phylogenetic Diversity Index
xdata = read.delim("metadata_alphadiv_socint_para_feed.txt", header = T, sep = "\t")

###### Check distributions for all the covariates ######
#Create object with covariates of interest
test.vars= c("ffr.prop", "fle.prop","rain", "ffl.prop","age_months", "soc.int", "richness.para")

#Create create a new vector to plot them all at once.
par(mfrow= c(3,3))  #Split the windows

#Loop for creating plots.
for(i in 1: length(test.vars)) {
  hist(xdata[,test.vars[i]], main =test.vars[i])}

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

#Social interactions
#Check for zeros
sum(xdata$soc.int==0)
#Compare distributions of original and log transform data
hist(xdata$soc.int)
hist(log(xdata$soc.int + 1))
#No changes after log transform.

#The distribution of feeding proportions, social interactions, parasite richness and precipitation is similar after log transform, so leave the 
#covariates as they are.

###Check number of data points per categorical variable -> look for outliers & number of observations.

#Data points per individual
plot(table(xdata$individual))
#Data points per sex
plot(table(xdata$sex))
#Data points per age_category
plot(table(xdata$age_category))
#Few data points for juveniles & infants
plot(table(xdata$group))

#######Correlation#####
# Check for correlations between covariates through scatterplots
plot(xdata$fle.prop, xdata$ffr.prop)
plot(xdata$rain, xdata$fle.prop)
plot(xdata$rain, xdata$ffr.prop)
plot(xdata$fle.prop, xdata$ffl.prop)
plot(xdata$rain, xdata$ffl.prop)
plot(xdata$soc.int, xdata$ffr.prop)
plot(xdata$soc.int, xdata$fle.prop)
plot(xdata$soc.int, xdata$ffl.prop)
plot(xdata$soc.int, xdata$rain)
plot(xdata$richness.para, xdata$fle.prop)
plot(xdata$richness.para, xdata$ffl.prop)
plot(xdata$richness.para, xdata$ffr.prop)
plot(xdata$richness.para, xdata$rain)

#Check for correlations between covariates using Pearson correlation
test.vars = c("ffr.prop", "fle.prop","rain", "ffl.prop", "age_months", "soc.int", "richness.para")
res=cor(xdata[,test.vars], method = "pearson", use = "complete.obs")
#Get the largest absolute correlation
max(abs(res[upper.tri(res)]))
#round to three decimals.
round(res, 3)
#               ffr.prop fle.prop   rain ffl.prop age_months soc.int richness.para
#ffr.prop         1.000   -0.141  0.010   -0.080     -0.017  -0.173        -0.103
#fle.prop        -0.141    1.000 -0.040    0.013     -0.098  -0.042         0.084
#rain             0.010   -0.040  1.000    0.041      0.033  -0.078        -0.082
#ffl.prop        -0.080    0.013  0.041    1.000     -0.015  -0.069         0.030
#age_months      -0.017   -0.098  0.033   -0.015      1.000  -0.154         0.104
#soc.int         -0.173   -0.042 -0.078   -0.069     -0.154   1.000        -0.171
#richness.para   -0.103    0.084 -0.082    0.030      0.104  -0.171         1.000

#No variables correlate as the highest value obtained is 0.1731673.

##### Random slopes ####
#Include random effect of individual and SampleID
#Test random slopes for feeding proportions, social interactions, parasite richness, age and precipitation

#Random slopes were determined by using a function developed by Dr. Roger Mundry.

source("diagnostic_fcns.r")

#Random slopes should be considered differently for factors & covariates.
#For a factor: include those with 2  or more observations per level of random effect.
#For a covariate: include those 3 or more unique observations of the fixed effect per level of the random effect. 
#IF the covariate does not meet the rule as a covariate, check if it meets it for a factor.
#If only some levels meet the requirement Roger recommends to include the random slope when at least half of the levels from the random effect do fulfill the rules.

#Example for interpretations
#log.affr_rate_within_individual (covariate)`
#9  10  11  12 tot -> number of unique observations (here should be 3 or more)
#1   1   2   1   5 -> levels of the random effect

#log.affr_rate_within_individual (covariate as factor)`
#0   1 tot -> number of unique observations (here should be 2 or more)
#2   3   5 -> levels of the random effect

#Model for checking random slopes using function fe.re.tab
xx.fe.re = fe.re.tab(
  fe.model = "alpha.div ~ sex + group + log.age_months + soc.int + ffr.prop + fle.prop + rain + ffl.prop + richness.para", #include her covariates to check for random slopes
  re = "(1|individual)", #include here the random effects
  data = xdata, other.vars = c("age_months")) #test age_months before log transform

xx.fe.re$summary

#According to the results the random slopes of feeding proportions, social interactions, parasite richness, age, and precipitationin individual should be included.
#Sex and group do not have random slopes.

#Extract the mean and standard deviation for all predictors before z transformation.
source("helpers.r")
t.data = xx.fe.re$data
wt.txt(c.descr(t.data$log.age_months))
wt.txt(c.descr(t.data$soc.int))
wt.txt(c.descr(t.data$ffr.prop))
wt.txt(c.descr(t.data$fle.prop))
wt.txt(c.descr(t.data$rain))
wt.txt(c.descr(t.data$ffl.prop))
wt.txt(c.descr(t.data$richness.para))

#z transform the covariates to ease model converge.
t.data$z.fle.prop=as.vector(scale(t.data$fle.prop))
t.data$z.rain = as.vector(scale(t.data$rain))
t.data$z.ffr.prop = as.vector(scale(t.data$ffr.prop))
t.data$z.ffl.prop = as.vector(scale(t.data$ffl.prop))
t.data$z.log.age_months = as.vector(scale(t.data$log.age_months))
t.data$z.soc.int = as.vector(scale(t.data$soc.int))
t.data$z.richness.para = as.vector(scale(t.data$richness.para))
str(t.data)

#Check distribution of the response variable alpha diversity (alpha.div)
hist(xdata$alpha.div)
#Log transform 
hist(log(xdata$alpha.div))
#After log transformation the distribution is not symmetrical.
#Transformation of alpha diversity ->Box-Cox power transformation.
xx=lm(alpha.div ~ sex + group + z.log.age_months + z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop + z.soc.int + z.richness.para,
      data=t.data)

#Box-Cox transformation####
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
full.wc = lmer(tr.alpha.div ~ sex + group + z.log.age_months + z.soc.int + z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop + z.richness.para+
                 (1 + z.log.age_months + z.soc.int + z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop + z.richness.para|individual), REML=F,
               data=t.data, control = contr)

diagnostics.plot(full.wc)

#Residuals against fitted values are not perfect but still acceptable.

#Revise results from the model
summary(full.wc)$varcor

#Model without correlations #####

full.woc = lmer(tr.alpha.div ~ sex + group + z.log.age_months + z.soc.int + z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop + z.richness.para+
                  (1 + z.log.age_months + z.soc.int + z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop + z.richness.para||individual), REML=F,
                data=t.data, control = contr)

#Compare log likelihoods of both models to determine if correlations should be kept.
logLik(full.wc)# -2542.677 (df=49)
logLik(full.woc) #-2548.906 (df=21)
#Keep correlations as there is a big difference between both models.

#Check for collinearity between predictors, this happens when values are higher than 5.
vif(full.wc)
vif(full.wc)[,3]^2
#No collinearity as the highest value was 1.439125 for parasite richness.

#Check distribution of random effects -> BLUPs (Best linear unbiased predictors).
ranef.diagn.plot(full.wc)
#They are acceptable.

#### Model stability ####
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
#They are acceptable with the exception of the individual and feeding on flowers.

### Inference ####
summary(full.wc)

##### Full-null model comparison #####

###Null model ####
#Prepare a null model excluding only the predictors of interest group, social interactions, age, sex and parasite richness from 
#the fixed effects but including their random slopes.

#Change optimizer to ease model converge
contr=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000))

null = lmer(tr.alpha.div ~ z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop +
              (1 + z.log.age_months + z.soc.int + z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop + z.richness.para|individual), REML=F,
            data=t.data, control = contr)

#Comparison between full-null model.
result = as.data.frame(anova(null, full.wc, test="Chisq"))
#Results shown in supplementary table S8
#p value  = 0.08830047
#The full and the null models are not significantly different.
#Export results

write.table(result, file ="full.null.alldataset.csv", sep = ",", col.names = T, row.names = T)

# Test the impact of single predictors ####

#Load package
library(lmerTest) #-> how to test single predictors in Gaussian models.

#Set optimizer to ease model convergence
contr.reml=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000))

#Calculate a new model in lmerTest for testing the effect of single predictors.
full.wc.reml = lmer(tr.alpha.div ~ sex + group + z.log.age_months + z.soc.int + z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop + z.richness.para+
                      (1 + z.log.age_months + z.soc.int + z.ffr.prop + z.fle.prop + z.rain + z.ffl.prop + z.richness.para|individual), REML=F,
                    data=t.data, control = contr.reml)

#Estimates model calculated with lmerTest
round(summary(full.wc.reml)$coefficients,3)
#Results:
#                   Estimate Std. Error      df t value Pr(>|t|)
#(Intercept)        26.124      1.641  29.235  15.915    0.000
#sexmale             1.153      0.924  94.322   1.248    0.215
#groupB             -0.622      1.859  31.946  -0.335    0.740
#groupF             -5.333      1.957  46.097  -2.725    0.009*
#groupJ             -2.116      1.808  32.113  -1.170    0.251
#z.log.age_months   -0.055      0.466 178.530  -0.118    0.906
#z.soc.int           0.869      0.486  11.289   1.786    0.101
#z.ffr.prop         -0.614      0.523  29.117  -1.173    0.250
#z.fle.prop          2.243      0.412  59.801   5.447    0.000*
#z.rain             -0.985      0.546   6.504  -1.805    0.117
#z.ffl.prop         -0.169      0.637  10.635  -0.265    0.796
#z.richness.para    -1.023      0.649  23.098  -1.577    0.128

detach(name = "package:lmerTest") #important when lmer will be used again.

#Export results#
estimates = as.data.frame(round(summary(full.wc.reml)$coefficients,3))
#Export table for results shown in supplementary table S8
write.table(estimates, file = "estimates_alpha_para_fulldataset.csv", sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)

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
plot.full.wc = lmer(tr.alpha.div ~ sex.male + group.B + group.F + group.J + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.richness.para +
                      (1 + z.log.age_months + z.soc.int + z.rain + z.fle.prop + z.ffr.prop + z.ffl.prop + z.richness.para|individual), REML=F,
                    data=t.data, control = contr)


#Plot for feeding on leaves ##

#Extract confidence intervals for feeding on leaves
plot.boot.full = boot.lmer(m = plot.full.wc, discard.warnings = F, nboots = 1000, para = T, n.cores = 5,
                           resol = 1000, level = 0.95, use = list("z.fle.prop")) #change list () according to the predictors of interest
#Check intervals                                                                      
plot.boot.full$ci.fitted


##Plot group -> figure 2a.
#Set axis
range(t.data$alpha.div)
yaxt.labels=seq(20,60, by=10)
yaxt.at=(yaxt.labels^transform.para$lambda-1)/transform.para$lambda
yaxt.at=yaxt.at/max((t.data$alpha.div^transform.para$lambda-1)/transform.para$lambda)
yaxt.at=yaxt.at*max(t.data$alpha.div)

#Set confidence intervals
ci=boot.full$ci.fitted
colnames(ci)[colnames(ci)%in%c("lower.cl","upper.cl")]=c("lwr", "upr")

#Plots
xx = factor.int.plot (
  plot.data = t.data, factors = "group", coefs=fixef(full.wc), response = t.data$tr.alpha.div, link=c("identity"),
  conf.int=ci, yaxt.labels=yaxt.labels, yaxt.at=yaxt.at, ylab = "Alpha diversity", ylim=NULL, size.fac=1.5,
  factor.seqs=NULL, factor.labels=NULL, to.show=c("bubbles", "boxes"), weights=NULL, which.q=c(3),
  pch=19, pt.col=grey(level=0.5, alpha = 0.5), rect.col=NA, border.col=NULL, est.ci.col=par("fg"), est.ci.lty=1, quiet=F, average.response=F, log.y=F, 
  bg="white", median.col="black", median.lwd=2, fitted.lwd=2, quant.col=NULL,
  cex.lab=1, cex.lab.x=1, cex.axis=1, reset.par=T, percentile.col=NULL, percentile.lwd=NULL, my.par=NULL, cex.axis.lab=1, xlim=NULL,
  add.range=F, range.pch=4)


##Plot feeding on leaves -> figure 2b
#Set confidence intervals
ci=plot.boot.full$ci.fitted
#Set axis
range(t.data$fle.prop)
x.labels=c(0, 0.1, 0.2)
x.at=(x.labels-mean(t.data$fle.prop))/sd(t.data$fle.prop)
#Plot
no.int.plot (
  plot.data = t.data, covariate = "z.fle.prop", coefs=fixef(plot.full.wc), response = t.data$tr.alpha.div, link=c("identity"),
  grid.resol=NA, 
  xlab="Feeding on leaves", ylab="Alpha diversity", x.at=x.at, x.labels=x.labels, y.at=yaxt.at, y.labels=yaxt.labels, xlim=NULL, ylim=NULL,
  cex.lab=1, cex.axis=1, size.fac=1.5, col=grey(level=0.5, alpha = 0.5), pch=19, conf.int=ci, ci.type=c("lines", "area"), lty=2, lwd=1, pt.lwd=1,
  my.par=list(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15, las=1), quiet=T, reset.par=F, mean.fun=mean,
  weights=NULL)