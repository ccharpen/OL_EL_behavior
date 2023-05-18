library(careless)
library(stringr)
library(dplyr)

rm(list = ls())

# set working directory to current source file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# helper functions --------------------------------------------------------
zscore = function(x) {
  return((x - mean(x, na.rm = T))/sd(x, na.rm = T))
}

# extract questionnaire data with questionnaire name (complete cases only)
get_complete <- function(data, grep_str, subscale='', invertGrep=F) {
  if (invertGrep) {
    data = data[, ls(data)[!grepl(grep_str, ls(data))]]
  } else {
    data = data[, ls(data)[grepl(grep_str, ls(data))]]
  }
  data = data[complete.cases(data),]
  if (subscale!='') {
    data = data[, names(data)[grepl(subscale, scales[names(data)])]]
  }
  return(data)
}

# intra-individual response variability
get_surveyIRVs = function(CIE_df, all_data, questionnaires) {
  # compute for each questionnaire
  for (q in questionnaires) {
    q_data = get_complete(all_data, q)
    q = gsub("\\", "", q, fixed=TRUE)
    CIE_df[, paste0('survey.irv.', q)] = irv(q_data)[rownames(CIE_df)] # compute IRVs for current questionnaire
    CIE_df[, paste0('survey.irv.', q, '.z')] = zscore(CIE_df[, paste0('survey.irv.', q)]) # z-score IRVs for current questionnaire
  }
  # get min, max, and average IRV
  max.irv =  apply(CIE_df[, grepl('.*\\irv.*\\.z', names(CIE_df))], 1, function(x) max(x, na.rm = T))
  min.irv = apply(CIE_df[, grepl('.*\\irv.*\\.z', names(CIE_df))], 1, function(x) min(x, na.rm = T))
  avg.irv = apply(CIE_df[, grepl('.*\\irv.*\\.z', names(CIE_df))], 1, function(x) mean(x, na.rm = T))
  CIE_df[, paste0('survey.irv.', 'max', '.z')] = max.irv
  CIE_df[, paste0('survey.irv.', 'min', '.z')] = min.irv
  CIE_df[, paste0('survey.irv.', 'avg', '.z')] = avg.irv
  return(CIE_df)
}

# splithalf function adapated from careless package's even-odd index
splithalf <- function(x, factors, diag = FALSE) {
  #initialize a result dataset
  #warning("Computation of even-odd has changed for consistency of interpretation
  #        with other indices. This change occurred in version 1.2.0. A higher
  #        score now indicates a greater likelihood of careless responding. If
  #        you have previously written code to cut score based on the output of
  #        this function, you should revise that code accordingly.")
  if(length(factors) == 1) {
    stop("You have called even-odd with only a single factor. \n The even-odd method requires multiple factors to work correctly.",
         call. = FALSE) }
  if(sum(factors) > ncol(x)) {
    stop("The number of items specified by 'factors' exceeds the number of columns in 'x'.",
         call. = FALSE) }
  if(sum(factors) != ncol(x)) {
    warning("The number of items specified by 'factors' does not match the number of columns in 'x'. \n Please check if this is what you want.",
            call. = FALSE) }
  
  # initalize empty list for persons holding the persons even scores and odd scores
  eo_vals <- vector("list", nrow(x))
  
  # Loop through each Person
  for(i in 1:nrow(x)) {
    # Initialize an object to hold the factor e/o means for the current person
    f <- matrix(rep(NA, 2*length(factors)), length(factors), ncol=2)
    start <- 1
    
    # loop through each factor
    for(j in 1:length(factors)) {
      if(j>1) start <- start + (factors[j-1])
      end <- (factors[j]-1) + start
      
      # Subset x with items for the current factor
      s <- x[i,start:end]
      ind <- seq(1:length(colnames(s)))
      e_ind <- which(ind %% 2 == 0)
      o_ind <- which(ind %% 2 == 1)
      f[j,1] <- mean(t(s[e_ind]), na.rm = TRUE)
      f[j,2] <- mean(t(s[o_ind]), na.rm = TRUE)
    }
    # assign the even and odd values to eo_vals
    eo_vals[[i]] <- f
  }
  
  #calculate number of even/odd pairs for which no comparison can be computed because of NAs
  eo_missing <- lapply(eo_vals, function(i) sum(!is.na(apply(i, 1, sum))))
  
  # scan for persons for which no even-odd can be calculated when all values are same, leading to
  # a correlation of NA because there is no variance/standard deviation.
  eo_sdzero <-  lapply(eo_vals, function(i) apply(i, 2, stats::sd))
  eo_sdzero <- sapply(eo_sdzero, function(i) any(i == 0, na.rm = TRUE))
  if(any(eo_sdzero, na.rm = TRUE)) warning("One or more observations have zero variance in even and/or odd values. \nThis results in NA values for these observations.\nIncluding more factors may alleviate this issue.",
                                           call. = FALSE)
  # Calculate within-person correlation between even and odd sub-scales
  # then apply the Spearman-Brown correction for split-half reliability
  # and store the result in the output vector.
  eo_cor <- sapply(eo_vals, function(f) {
    # suppressWarnings for standard deviation 0 which happens when each value-pairs is same
    val <- suppressWarnings(stats::cor(f[, 1], f[, 2], use = "pairwise.complete.obs"))
    val <- (2 * val) / (1 + val) #split-half
    if (!is.na(val) && val < -1) val <- -1
    return(val)
  })
  
  # transform eo  such that higher scores are indicating carelessness
  # eo_cor = 0 - eo_cor
  if(diag == FALSE) {return(eo_cor)}
  else {return(data.frame(eo_cor, eo_missing))}
}

# repeatedly call splithalf, computing rolling average
# we may want to edit this so that we get a split half reliability score for each questionnaires?
repeat_splithalf <- function(data, scales_list, name='', niter=100) {
  #res_fname = './splithalf_average.RDS'
  lengths = unlist(lapply(scales_list,length))
  for (it in 1:niter) {
    scales_list = lapply(scales_list,sample)
    scales = unlist(scales_list)
    if (name=='') {
      sh = splithalf(data[, scales], factors=lengths)
    } else{
      sh = splithalf(data[, paste(name, scales, sep='_')], factors=lengths)
    }
    if (it==1) {
      results = sh
    } else {
      results = results + (sh-results)/it
    }
    res = as.data.frame(results, row.names=rownames(data))
    #saveRDS(res,res_fname)
  }
  return(res)
}

BDI = c(1:21)
State = c(1:20)
Trait = c(1:20)
SRS = c(1:65)
LSAS = c(1:48)
scale_items_dict = list(BDI=BDI, State=State, Trait=Trait, SRS=SRS, LSAS=LSAS)

# get splithalf reliability (correlation of scores when splitting each scale/subscale into random halves)
get_surveySplitHalfRels = function(CIE_df, all_data, questionnaires=c()) {
  questionnaires_all = unlist(unique(lapply(str_split(names(all_data),'_'), function(x) x[1])))
  questionnaires_new = questionnaires_all[questionnaires_all %in% gsub("\\", "", questionnaires,fixed=T)]
  if (length(questionnaires) != length(questionnaires_new)) {
    print(questionnaires)
  }
  questionnaires = questionnaires_new
  scale_items_all = list()
  qI = 1
  sI = 1
  q_added = list(names=list(), sI=list(), qI=list())
  # build up scale_items_all which contains information about which questionnaires belong to which scale for use by splithalf function
  for (qname in questionnaires) {
    q_data = get_complete(all_data, qname)
    if (qname %in% names(scale_items_dict)) {
      scale_items_list = scale_items_dict[[qname]]
      
      q_added$names = append(q_added$names, qname)
      q_added$qI = append(q_added$qI, qI)
      q_added$sI = append(q_added$sI, sI)
      scale_items_all[[sI]] = qI+scale_items_list-1
      sI = sI + 1
      qI = qI + length(unlist(scale_items_list))
    } else {
      scale_items_all[[sI]] = seq(qI, qI+ncol(q_data)-1)
      qI = qI + ncol(q_data)
      sI = sI + 1
    }
  }
  all_data = all_data[, unlist(lapply(str_split(names(all_data),'_'), function(x) x[1])) %in% questionnaires]
  sh_res = repeat_splithalf(all_data, scale_items_all)
  CIE_df[, 'survey.shr'] = sh_res[rownames(CIE_df),]
  CIE_df[, 'survey.shr.z'] = zscore(CIE_df[, 'survey.shr'])
  return(CIE_df)
}

# get correlation of pyschometric synonyms  (r > 0.7 are considered synonyms)
get_surveyPsychSyn = function(CIE_df, all_data, questionnaires) {
  all_data = all_data[,grepl(paste(questionnaires, collapse='|'), names(all_data))]
  all_syn <- data.frame(psychsyn(all_data, critval = 0.65, diag = T), row.names=rownames(all_data))
  CIE_df[rownames(all_syn), 'survey.syn.r'] = all_syn$cor
  CIE_df[, 'survey.syn.r.z'] = zscore(CIE_df[, 'survey.syn.r'])
  return(CIE_df)
}

# get correlation of pyschometric antonyms (r < -0.45 are considered antonyms)
get_surveyPsychAnt = function(CIE_df, all_data, questionnaires) {
  all_data = all_data[,grepl(paste(questionnaires, collapse='|'), names(all_data))]
  all_ant <- data.frame(psychant(all_data, critval = 0.5, diag = T), row.names=rownames(all_data))
  CIE_df[rownames(all_ant), 'survey.ant.r'] = all_ant$cor
  CIE_df[, 'survey.ant.r.z'] = zscore(CIE_df[, 'survey.ant.r'])
  return(CIE_df)
}


# read data --------------------------------------------------------------------------------

csv_fname <- '../../data/total_questionnaire_items.csv'      # csv file to hold questionnaire data (reversed coded item scores)
all_data = read.csv(csv_fname)
rownames(all_data) <- all_data$uniqueID
all_data$uniqueID <- NULL
all_data$subNb <- NULL
all_data$study <- NULL

csv_raw_fname <- '../../data/total_questionnare_items_raw.csv'      # csv file to hold questionnaire data (raw item scores)
all_data_raw = read.csv(csv_raw_fname)
rownames(all_data_raw) <- all_data_raw$uniqueID
all_data_raw$uniqueID <- NULL
all_data_raw$subNb <- NULL
all_data_raw$study <- NULL

# get C/IE measures -------------------------------------------------------
CIE_df = data.frame(row.names = rownames(all_data))
CIE_raw_df = data.frame(row.names = rownames(all_data_raw))

# questionnaires of interest
questionnaire_names = c('BDI', 'State', 'Trait', 'SRS', 'LSAS')

# get intra-individual response variability
CIE_df = get_surveyIRVs(CIE_df, all_data, questionnaire_names)
CIE_raw_df = get_surveyIRVs(CIE_raw_df, all_data_raw, questionnaire_names)

# get split-half reliability (100 random partitions of dataset by default)
CIE_df = get_surveySplitHalfRels(CIE_df, all_data, questionnaire_names)
CIE_raw_df = get_surveySplitHalfRels(CIE_raw_df, all_data_raw, questionnaire_names)


# histogram for psychometric synonyms and antonyms
hist(psychsyn_critval(all_data_raw)$Freq,100,
     main = expression('Histogram of' ~ italic('r') ~ 'values'),
     xlim = c(-1,1))
abline(v = c(-0.5, 0.65), lty='dashed', col='red')
length(which(psychsyn_critval(all_data_raw)$Freq>(0.65)))
length(which(psychsyn_critval(all_data_raw)$Freq<(-0.5)))

# get psychometric synonyms 
CIE_raw_df = get_surveyPsychSyn(CIE_raw_df, all_data_raw, questionnaire_names)

# get psychometric and antonyms
CIE_raw_df = get_surveyPsychAnt(CIE_raw_df, all_data_raw, questionnaire_names)

df_total = data.frame(IRV_STAI_State = CIE_raw_df$survey.irv.State,
                      IRV_STAI_State_z = CIE_raw_df$survey.irv.State.z,
                      IRV_STAI_Trait = CIE_raw_df$survey.irv.Trait,
                      IRV_STAI_Trait_z = CIE_raw_df$survey.irv.Trait.z,
                      IRV_SRS = CIE_raw_df$survey.irv.SRS,
                      IRV_SRS_z = CIE_raw_df$survey.irv.SRS.z,
                      SHR = CIE_df$survey.shr,
                      SHR_z = CIE_df$survey.shr.z,
                      Syn = CIE_raw_df$survey.syn.r,
                      Syn_z = CIE_raw_df$survey.syn.r.z,
                      Ant = CIE_raw_df$survey.ant.r,
                      Ant_z = CIE_raw_df$survey.ant.r.z,
                      BDI_Q1 = all_data$BDI_Q1,
                      LSAS_Q1 = all_data$LSAS_Q1_A,
                      row.names = rownames(all_data))

df_total$IRV_avg = rowMeans(subset(df_total, select = c(IRV_STAI_State_z, IRV_STAI_Trait_z, IRV_SRS_z)),na.rm=TRUE)

temp = read.csv('../../data/summary_data_pooled.csv')
temp = temp[temp$excluded_miss==0,]
rownames(temp) <- (1:length(temp$subNb))
df_catch = temp[c("study","subNb","catchQ_BDI","catchQ_SRS","catchQ_STAI")]

df_total = merge(df_total, df_catch, by='row.names')
names(df_total)[names(df_total) == 'Row.names'] <- 'uniqueID'
df_total$uniqueID <- as.numeric(df_total$uniqueID)
df_total <- df_total[order(df_total$uniqueID),]
df_total <- df_total %>% relocate(study, subNb)
rownames(df_total)<-NULL


df_total$catchQ_Total = rowSums(subset(df_total, select=c(catchQ_BDI, catchQ_SRS, catchQ_STAI)),na.rm=TRUE)
df_total$catchQ_Acc = rowMeans(subset(df_total, select=c(catchQ_BDI, catchQ_SRS, catchQ_STAI)),na.rm=TRUE)

df_total$excluded_catchQs = df_total$catchQ_Acc < 1
df_total$excluded_IRV_STAI_State = df_total$IRV_STAI_State_z < -2
df_total$excluded_IRV_STAI_Trait = df_total$IRV_STAI_Trait_z < -2
df_total$excluded_IRV_SRS = df_total$IRV_SRS_z < -2
df_total$excluded_IRV_avg = df_total$IRV_avg < -2
df_total$excluded_SHR = df_total$SHR_z < -2
df_total$excluded_Syn = df_total$Syn_z < -2
df_total$excluded_Ant = df_total$Ant_z > 2
df_total$excluded_na = (is.na(df_total$catchQ_Total) | is.na(df_total$IRV_STAI_State_z) | is.na(df_total$IRV_STAI_Trait_z) 
                        | is.na(df_total$IRV_SRS_z) | is.na(df_total$SHR_z) | is.na(df_total$Syn_z) | is.na(df_total$Ant_z)
                        | is.na(df_total$BDI_Q1) | is.na(df_total$LSAS_Q1))
df_total$excluded = (df_total$excluded_IRV_STAI_State | df_total$excluded_IRV_STAI_Trait
                     | df_total$excluded_IRV_SRS | df_total$excluded_SHR
                     | df_total$excluded_Syn | df_total$excluded_Ant | df_total$excluded_na)
df_total$excluded2 = (df_total$excluded_IRV_avg | df_total$excluded_SHR
                      | df_total$excluded_Syn | df_total$excluded_Ant | df_total$excluded_na)

# write output as csv   
write.csv(x=df_total, file='Careless_exclusions_output.csv')
