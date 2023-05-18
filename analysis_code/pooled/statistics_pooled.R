library(tidyverse)
library(rstatix)
library(lme4)
library(lmerTest)
library(nloptr)
library(emmeans)


rm(list = ls())

# set working directory to current source file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read data
df <- read.csv("recap_factor_groups.csv")

#set control parameters for mixed effects models
control_params = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))

#create a long format table for lmer
df_f <- df[c(1:9,83:90)]
#df_f <- na.omit(df_f)
dfl <- pivot_longer(df_f, cols = 10:17,  values_to = "score", names_to = "factor")

#convert group and study variables to factor
dfl[,'group'] <- lapply(dfl[,'group'],factor)
dfl[,'study'] <- lapply(dfl[,'study'],factor)

#run models
mod_group <- lmer(score ~ factor*group + study + (1|uniqueID), data = dfl, REML=FALSE, control = control_params)
anova(mod_group)

mod_group_cov <- lmer(score ~ factor*group + gender + age + education + ICAR_score + study + (1|uniqueID), data = dfl, REML=FALSE, control = control_params)
anova(mod_group_cov)
confint(mod_group_cov)

#extract marginal means
emm1 <- emmeans(mod_group_cov, specs = pairwise ~ factor | group, pbkrtest.limit = 4512)
emm1$emmeans
emm1$contrasts
eff_size(emm1, sigma = sigma(mod_group_cov), edf = 3983)

emm2 <- emmeans(mod_group_cov, specs = pairwise ~ group | factor, pbkrtest.limit = 4512)
emm2$emmeans
emm2$contrasts
eff_size(emm2, sigma = sigma(mod_group_cov), edf = 3983)

emm1_na <- emmeans(mod_group_cov, specs = pairwise ~ factor | group, pbkrtest.limit = 4512, adjust = 'none')
emm1_na$emmeans
emm1_na$contrasts

cor.test(df$F3_autism,df$F7_traitAnx)
cor.test(df[df$group==1,'F3_autism'],df[df$group==1,'F7_traitAnx'])
cor.test(df[df$group==2,'F3_autism'],df[df$group==2,'F7_traitAnx'])
cor.test(df[df$group==3,'F3_autism'],df[df$group==3,'F7_traitAnx'])
cor.test(df[df$group==4,'F3_autism'],df[df$group==4,'F7_traitAnx'])
cor.test(df[df$group==5,'F3_autism'],df[df$group==5,'F7_traitAnx'])

