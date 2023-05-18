library(psych)
library(dplyr)
library(nFactors)
library(corrplot)
library(ggplot2)
library(stringr)
library(gridExtra)
library(paran)
library(EFA.dimensions)
library(reshape2)

rm(list = ls())

# set working directory to current source file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# define a function for plotting factor loadings
plot_factor <- function(df_loadings, x, y, title) {
  p_f <- ggplot(df_loadings) +
    geom_bar(aes(x = x, y = y, fill=Questionnaire), colour="black", size=0.1, stat="identity") +
    geom_hline(yintercept = 0) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.size = unit(.6, "lines"),
          axis.line.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x.bottom =element_blank(),
          panel.grid = element_blank()) +
    ggtitle(title) +
    ylab("Loading Strength") +
    scale_fill_brewer(palette="Pastel1")
  return(p_f)
}

#function for calculating p-values of correlation matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# load the data
csv_fname = '../../data/total_questionnaire_items.csv'
all_data = read.csv(csv_fname)

# exclusion
exclusion = read.csv('Careless_exclusions_output.csv')
exc = exclusion[c("uniqueID","excluded2")]
all_data = merge(all_data, exc, by="uniqueID")

data = all_data[all_data$excluded2==FALSE,]
id_list = data$uniqueID
rownames(data) <- data$uniqueID
data$uniqueID <- NULL
data$study <- NULL
data$subNb <- NULL
data$excluded2 <- NULL

# correlation matrix & eigenvalues for factor analysis
corr_mat = cor(data, use="pairwise")
scree(corr_mat)
eigens = eigen(corr_mat)

# initialize data frame for comparing models with different number of factors
df_mod_comp <- as.data.frame(matrix(nrow=20,ncol=13))
colnames(df_mod_comp) <- c("nFactors","rms","crms","TLI","RMSEA","Fit","FitOff","BIC","chi2","df","chi2-diff","df-diff","sig_threshold_001")
for (nf in 1:20) {
  print(nf)
  fa_model_wls = fa(data, nfactors = nf, fm='wls', rotate = 'oblimin')
  if (nf==1) {
    df_mod_comp[nf,] = c(nf, fa_model_wls$rms, fa_model_wls$crms, fa_model_wls$TLI, as.numeric(fa_model_wls$RMSEA[1]), fa_model_wls$fit, fa_model_wls$fit.off,
                         fa_model_wls$BIC, fa_model_wls$chi, fa_model_wls$dof, NA, NA, NA)
  } else {
    df_mod_comp[nf,] = c(nf, fa_model_wls$rms, fa_model_wls$crms, fa_model_wls$TLI, as.numeric(fa_model_wls$RMSEA[1]), fa_model_wls$fit, fa_model_wls$fit.off,
                         fa_model_wls$BIC, fa_model_wls$chi, fa_model_wls$dof, chi_prev-fa_model_wls$chi, df_prev-fa_model_wls$dof, qchisq(0.001,df_prev-fa_model_wls$dof,lower.tail=FALSE))
  }
  chi_prev = fa_model_wls$chi
  df_prev = fa_model_wls$dof
}
write.csv(x=df_mod_comp, file='FA_model_comparison_wls.csv')
#TLI: Tucker Lewis fit index, typically reported in SEM. Generally want > .9
#RMSEA: Root mean square error of approximation. Also reported is the so-called 'test of close fit'
#FitOff: Fit based upon off diagonal values: you can think of it as 1 - resid^2 / cor^2, or a kind of R2 applied to a correlation matrix instead of raw data
#BIC: Useful for model comparison purposes only.
#SABIC: sample-size corrected BIC

#plot BIC
ggplot(df_mod_comp, mapping = aes(x=nFactors, y=BIC)) + geom_point() + geom_line()


### Run factor analysis with 8 factors
nf=8
fa_model = fa(data, nfactors = nf, fm='wls', rotate = 'oblimin')

df_loadings = data.frame(fa_model$loadings[,1:nf])
df_loadings <- df_loadings[,order(names(df_loadings))]
df_loadings$ q = 1:174
df_loadings$Questionnaire = factor(c(rep("BDI",21),rep("STAI-S",20),rep("STAI-T",20),rep("SRS",65),rep("LSAS",48)))
write.csv(x=df_loadings, file='factor_loadings_8F.csv')

# plot factor loadings and save figure
p_f1 <- plot_factor(df_loadings, df_loadings$q, df_loadings$WLS1, "Factor 1")
p_f2 <- plot_factor(df_loadings, df_loadings$q, df_loadings$WLS2, "Factor 2")
p_f3 <- plot_factor(df_loadings, df_loadings$q, df_loadings$WLS3, "Factor 3")
p_f4 <- plot_factor(df_loadings, df_loadings$q, df_loadings$WLS4, "Factor 4")
p_f5 <- plot_factor(df_loadings, df_loadings$q, df_loadings$WLS5, "Factor 5")
p_f6 <- plot_factor(df_loadings, df_loadings$q, df_loadings$WLS6, "Factor 6")
p_f7 <- plot_factor(df_loadings, df_loadings$q, df_loadings$WLS7, "Factor 7")
p_f8 <- plot_factor(df_loadings, df_loadings$q, df_loadings$WLS8, "Factor 8")

png(filename = 'loading_strength_excl_8F.png', width = 2000, height = 1200, res = 200)
grid.arrange(p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, p_f8, ncol=2)
dev.off()

# extract individual subject factor scores
indiv_scores = data.frame(fa_model$scores)
indiv_scores <- indiv_scores[,order(names(indiv_scores))]
indiv_scores$uniqueID = id_list
write.csv(x=indiv_scores, file='factor_scores_excl_8F.csv')

#correlation matrix of loadings across items
M <- cor(df_loadings[,1:nf])
p.mat <- cor.mtest(df_loadings[,1:nf])
corrplot(M, method="color", type="upper", addCoef.col = "black", tl.col="black", tl.srt=45,
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", diag=FALSE)

#correlation matrix of scores across subjects
M <- cor(indiv_scores[,1:nf])
p.mat <- cor.mtest(indiv_scores[,1:nf])
corrplot(M, method="color", type="upper", addCoef.col = "black", tl.col="black", tl.srt=45,
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", diag=FALSE)


