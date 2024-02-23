library(tidyverse)
library(rstatix)
library(lme4)
library(lmerTest)
library(nloptr)
library(ggplot2)

rm(list = ls())

#set working directory to current source file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read data
df <- read.csv("recap_analyses_study2.csv")

#set control parameters for mixed effects models
control_params = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))


## OL CHOICE
#create a long format table for lmer
df_OLch <- df[c(1:6,16:23)]
df_OLch <- na.omit(df_OLch)
dfl <- pivot_longer(df_OLch, cols = 7:14,  values_to = "OLch", names_to = c("OLunc","ELunc","Mag"), names_pattern = "OLch_O(.)uE(.)uM(.)")

#convert group variable to factor
dfl[,'group'] <- lapply(dfl[,'group'],factor)

#mixed-effect models, all subjects (Fig S2B)
mod <- lmer(OLch ~ OLunc*ELunc*Mag + (1|subNb), data = dfl, REML=FALSE, control = control_params)
anova(mod)

#mixed-effect models, interaction with group, controlling for covariates
mod_g_cov <- lmer(OLch ~ OLunc*ELunc*Mag*group + gender + age + education + ICAR_score + (1|subNb), data = dfl, REML=FALSE, control = control_params)
anova(mod_g_cov)

#repeat averaging across magnitude (Fig 3B)
mod <- lmer(OLch ~ OLunc*ELunc + (1|subNb), data = dfl, REML=FALSE, control = control_params)
anova(mod)
confint(mod)



## LEARNING CURVES
df_learn <- df[c(1:6,24:31)]
df_learn <- na.omit(df_learn)
dfl6 <- pivot_longer(df_learn, cols = 7:14,  values_to = "accuracy", names_to = c("trial"), names_pattern = "LearningT(.)")
dfl6[,'group'] <- lapply(dfl6[,'group'],factor)
dfl6 <- transform(dfl6, trial = as.numeric(trial))

#Fig 2B
mod_nocov <- lmer(accuracy ~ trial + (1|subNb), data = dfl6, REML=FALSE, control = control_params)
anova(mod_nocov)
confint(mod_nocov)

#Fig 5D
mod_group_cov <- lmer(accuracy ~ trial*group + gender + age + education + ICAR_score + (1|subNb), data = dfl6, REML=FALSE, control = control_params)
anova(mod_group_cov)
confint(mod_group_cov)



## GLME MAIN EFFECTS
df_main <- df[c(1:6,34:35)]
df_main <- na.omit(df_main)
dfl2 <- pivot_longer(df_main, cols = 7:8,  values_to = "glme_effect", names_to = c("strategy"), names_pattern = "GLME_(..)")
dfl2[,'group'] <- lapply(dfl2[,'group'],factor)

#Fig 2F
mod_nocov <- lmer(glme_effect ~ strategy + (1|subNb), data = dfl2, REML=FALSE, control = control_params)
anova(mod_nocov)
confint(mod_nocov)

#Fig 6C-D
mod_group_cov <- lmer(glme_effect ~ strategy*group + gender + age + education + ICAR_score + (1|subNb), data = dfl2, REML=FALSE, control = control_params)
anova(mod_group_cov)
confint(mod_group_cov)

t.test(df[df$group==2,"GLME_EL"],df[df$group==2,"GLME_OL"], paired = TRUE, var.equal = TRUE)
t.test(df[df$group==3,"GLME_EL"],df[df$group==3,"GLME_OL"], paired = TRUE, var.equal = TRUE)


## GLME effects - interaction with OL uncertainty
df_OLU <- df[c(1:6,36:39)]
df_OLU <- na.omit(df_OLU)
dfl3 <- pivot_longer(df_OLU, cols = 7:10,  values_to = "glme_effect", names_to = c("strategy","OLunc"), names_pattern = "GLME_(..)_([A-Za-z]+)OLU")
dfl3[,'group'] <- lapply(dfl3[,'group'],factor)

#Fig 3D
mod_nocov <- lmer(glme_effect ~ strategy*OLunc + (1|subNb), data = dfl3, REML=FALSE, control = control_params)
anova(mod_nocov)
confint(mod_nocov)

#Fig 7E
mod_group_cov <- lmer(glme_effect ~ strategy*OLunc*group + gender + age + education + ICAR_score + (1|subNb), data = dfl3, REML=FALSE, control = control_params)
anova(mod_group_cov)

diff_dyn = df[df$group==5,"GLME_OL_LowOLU"] - df[df$group==5,"GLME_OL_HighOLU"]
diff_fix = df[df$group==4,"GLME_OL_LowOLU"] - df[df$group==4,"GLME_OL_HighOLU"]
t.test(diff_dyn, diff_fix)

diff_dyn = df[df$group==5,"GLME_EL_LowOLU"] - df[df$group==5,"GLME_EL_HighOLU"]
diff_fix = df[df$group==4,"GLME_EL_LowOLU"] - df[df$group==4,"GLME_EL_HighOLU"]
t.test(diff_dyn, diff_fix)


## GLME effects - interaction with EL uncertainty
df_ELU <- df[c(1:6,40:43)]
df_ELU <- na.omit(df_ELU)
dfl4 <- pivot_longer(df_ELU, cols = 7:10,  values_to = "glme_effect", names_to = c("strategy","ELunc"), names_pattern = "GLME_(..)_([A-Za-z]+)ELU")
dfl4[,'group'] <- lapply(dfl4[,'group'],factor)

#Fig 3D
mod_nocov <- lmer(glme_effect ~ strategy*ELunc + (1|subNb), data = dfl4, REML=FALSE, control = control_params)
anova(mod_nocov)
confint(mod_nocov)

#Fig 7F
mod_group_cov <- lmer(glme_effect ~ strategy*ELunc*group + gender + age + education + ICAR_score + (1|subNb), data = dfl4, REML=FALSE, control = control_params)
anova(mod_group_cov)

diff_dyn = df[df$group==5,"GLME_OL_LowELU"] - df[df$group==5,"GLME_OL_HighELU"]
diff_fix = df[df$group==4,"GLME_OL_LowELU"] - df[df$group==4,"GLME_OL_HighELU"]
t.test(diff_dyn, diff_fix)

diff_dyn = df[df$group==5,"GLME_EL_LowELU"] - df[df$group==5,"GLME_EL_HighELU"]
diff_fix = df[df$group==4,"GLME_EL_LowELU"] - df[df$group==4,"GLME_EL_HighELU"]
t.test(diff_dyn, diff_fix)



## GLME effects - interaction with Magnitude
df_mag <- df[c(1:6,44:47)]
df_mag <- na.omit(df_mag)
dfl5 <- pivot_longer(df_mag, cols = 7:10,  values_to = "glme_effect", names_to = c("strategy","Mag"), names_pattern = "GLME_(..)_([A-Za-z]+)Mag")
dfl5[,'group'] <- lapply(dfl5[,'group'],factor)

#Fig S2D
mod_nocov <- lmer(glme_effect ~ strategy*Mag + (1|subNb), data = dfl5, REML=FALSE, control = control_params)
anova(mod_nocov)
confint(mod_nocov)

#Fig S2F
mod_group_cov <- lmer(glme_effect ~ strategy*Mag*group + gender + age + education + ICAR_score + (1|subNb), data = dfl5, REML=FALSE, control = control_params)
anova(mod_group_cov)



# GROUP DIFFERENCES IN 1 VARIABLE
anova(lm(baseIndex ~ factor(group), data = df))
anova(lm(arbIndex ~ factor(group), data = df)) #Fig 7B
anova(lm(abs(colorBias) ~ factor(group), data = df))
anova(lm(abs(handBias) ~ factor(group), data = df))
anova(lm(abs(stickyAct) ~ factor(group), data = df))
anova(lm(abs(actImit) ~ factor(group), data = df))
anova(lm(OL_alpha ~ factor(group), data = df))
anova(lm(EL_alpha ~ factor(group), data = df))
anova(lm(EL_magBoost ~ factor(group), data = df))
anova(lm(wOLEL_fix ~ factor(group), data = df))
anova(lm(biasOLEL_dyn ~ factor(group), data = df))
anova(lm(age ~ factor(group), data = df)) #Fig S8E
anova(lm(education ~ factor(group), data = df)) #Fig S8G
anova(lm(ICAR_score ~ factor(group), data = df)) #Fig S8H

#Fig S8F
chisq.test(table(df[!(df$gender=="N"),c("group","gender")])) #compare gender distribution ignoring non-binary cases
chisq.test(table(df[,c("group","gender")]))

#controlling for covariates
lm_arb <- lm(arbIndex ~ factor(group) + gender + age + education + ICAR_score, data = df)
anova(lm_arb)
confint(lm_arb)

t.test(df[df$group==5,"arbIndex"],df[df$group==4,"arbIndex"]) #Fig 7B
cohens_d(df[df$group>=4,], arbIndex ~ group)


## Arbitration weight (Fig S6C)
df_w <- df[c(1:6,62:69)]
#mixed-effect models, all subjects
dfl <- pivot_longer(df_w, cols = 7:14,  values_to = "w", names_to = c("OLunc","ELunc","Mag"), names_pattern = "w_O(.)uE(.)uM(.)")
mod <- lmer(w ~ OLunc*ELunc*Mag + (1|subNb), data = dfl, REML=FALSE, control = control_params)
anova(mod)

#mixed-effect models, interaction with group, controlling for covariates
#exclude missing values
df_w <- na.omit(df_w)
dfl <- pivot_longer(df_w, cols = 7:14,  values_to = "w", names_to = c("OLunc","ELunc","Mag"), names_pattern = "w_O(.)uE(.)uM(.)")
dfl[,'group'] <- lapply(dfl[,'group'],factor)
mod_g_cov <- lmer(w ~ OLunc*ELunc*Mag*group + gender + age + education + ICAR_score + (1|subNb), data = dfl, REML=FALSE, control = control_params)
anova(mod_g_cov)

