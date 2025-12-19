# OL_EL_behavior
Data and code for study investigating heterogeneity in arbitration behavior between experiential and observational learning, and its relevance for psychopathology.
Manuscript is published at:

Charpentier, C.J., Wu, Q., Min, S. et al. Heterogeneity in strategy use during arbitration between experiential and observational learning. Nat Commun 15, 4436 (2024). https://doi.org/10.1038/s41467-024-48548-y

## task_code folder
This folder contains the code to run the task online and uses jsPsych and HTML. A demo of the task instructions, practice, and first 10 trials is available at: https://obsexplearn.web.app/

## data folder
This folder contains the raw and summarized data for each of the two studies.
- The files 'data_study1.csv' and 'data_study2.csv' contain trial-level task data, including information about the trial sequence and design, choices, reaction times and outcomes.
- The file 'summary_data_pooled.csv' contains participant-level data, including study date/time, demographics, overall reaction times, performance metrics, bonus calculation and summary questionnaire scores.
- The files 'total_questionnaire_items.csv' and 'total_questionnaire_items_raw.csv' contain participant-level individual item scores from questionnaires, used in the careless responding analysis and factor analysis. The 'raw' item scores are the raw responses on the likert scale for each questionnaire (used for identifying patterns of careless responding), while the other file contains reverse-coded items (used for scoring).

## analysis_code folder
This folder contains all the code, as well as some output files, needed to run the analyses, statistics, and plot the findings. The analyses were mostly performed in MATLAB (R2020b), while some mixed-effect models and statistics were computed in R (version 4.1.2).
The 'model_recovery.m' is a script that perform recovery across the 5 computational models and plots a confusion matrix.

### dependencies sub-folder
This folder contains various functions needed to run the code:
- cbm-master: the toolbox (_cbm_, see Piray et al, 2019 and github repo: https://github.com/payampiray/cbm) used for fitting computational models
- model_functions: the functions used to compute the log-likelihood associared with each model, as well as generating choice data
- plotSpread: a toolbox for plotting individual data point distributions
- redwhiteblue.m: a function to plot correlation maps in a blue to white to red scale
- Violin.m and violinplot.m: functions to make violin plots (used in Figure 8)

### study1 & study2 sub-folders
These folders contain the different steps of the analysis, and corresponding output files, for each step:
- 'behavior_glm' performs behavioral and glm analyses of trial-by-trial data, and outputs the 'Behavioral_variables.mat' file
- 'model_fitting' performs model fitting, and parameter recovery, and outputs files contained in the 'model_fitting_outputs' folder
- 'group_analyses' performs group-level analyses, whereby behavioral metrics are computed and compared across the 5 groups identified during model fitting
- The R script 'statistics_final' performs various statistics and mixed-effects models reported in the paper
- The 'recap_analyses' csv file summarizes behavioral and model-derived variables at the participant level, used for plots and statistics

### pooled sub-folder
This folder performs individual difference analyses related to psychiatric symptom dimension measures, pooled across both studies given their exploratory nature.
- The R script 'careless_exclusions.R' performs analyses of careless responding in order to exclude participants who answered the self-report questionnaires carelessly. The output of these careless measures, and resulting exclusions, are reported in the 'Careless_exclusions_output.csv' file.
- The R script 'factor_analysis_final.R' performs factor analyses over the individual questionnaire items to extract relevant and separate psychiatric symptom dimensions. It outputs a file comparing factor analysis models with different number of factors ('FA_model_comparison_wls.csv'), the factor loadings for each item ('factor_loadings_8F.csv') and the individual participant's factor scores ('factor_scores_excl_8F.csv').
- Group-level analyses, i.e. comparing factor scores across the different groups, and the resulting plots, are performed in the 'group_differences_factors.m' script, with corresponding statistics run in R ('statistics_pooled.R')

## Reproducibility of Figures and Tables
Figures and tables presented in the paper (main text and supplementary information) can be reproduced from the following scripts contained in the analysis_code folder:
#### Main text
- Table 1: AIC, Frequency and Nbest (%best) columns are computed in /study1/model_fitting_study1.m and /study2/model_fitting_study2.m; OOS accuracy column is computed in /study1/out_of_sample_accuracy_study1.m and /study2/out_of_sample_accuracy_study2.m
- Figure 1B: /study1/behavior_glm_study1.m
- Figure 2: /study1/behavior_glm_study1.m and /study2/behavior_glm_study2.m
- Figure 3: /study1/behavior_glm_study1.m and /study2/behavior_glm_study2.m
- Figure 4: /study1/posterior_predictive_checks_study1.m and /study2/posterior_predictive_checks_study2.m
- Figure 5: /study1/out_of_sample_accuracy_study1.m and /study2/out_of_sample_accuracy_study2.m for Figure 5A-B; /study1/posterior_predictive_checks_study1.m and /study2/posterior_predictive_checks_study2.m for Figure 5C-D
- Figure 6: /study1/posterior_predictive_checks_study1.m and /study2/posterior_predictive_checks_study2.m
- Figure 7: /study1/group_analyses_study1.m and /study2/group_analyses_study2.m
- Figure 8: /pooled/group_differences_factors.m
#### Supplementary information
- Table S1: /study1/behavior_glm_study1.m and /study2/behavior_glm_study2.m
- Table S2: /study1/behavior_glm_study1.m and /study2/behavior_glm_study2.m
- Table S3: AIC column is computed in /study1/model_fitting_study1.m and /study2/model_fitting_study2.m; OOS accuracy column is computed in /study1/out_of_sample_accuracy_study1.m and /study2/out_of_sample_accuracy_study2.m
- Table S4: /study1/statistics_s1_final.R and /study2/statistics_s2_final.R
- Table S5: /pooled/statistics_pooled.R
- Figure S1: /study1/behavior_glm_study1.m and /study2/behavior_glm_study2.m
- Figure S2: /study1/behavior_glm_study1.m and /study2/behavior_glm_study2.m for Figure S2A-D; /study1/group_analyses_study1.m and /study2/group_analyses_study2.m for Figure S2E-F
- Figure S3: /model_recovery.m for Figure S3A; /study1/model_fitting_study1.m for Figure S3B-F
- Figure S4: /study1/posterior_predictive_checks_study1.m and /study2/posterior_predictive_checks_study2.m
- Figure S5: /study1/posterior_predictive_checks_study1.m and /study2/posterior_predictive_checks_study2.m
- Figure S6: /study1/model_fitting_study1.m and /study2/model_fitting_study2.m
- Figure S7: /study1/group_analyses_study1.m and /study2/group_analyses_study2.m
- Figure S8: /study1/posterior_predictive_checks_study1.m and /study2/posterior_predictive_checks_study2.m
- Figure S9: /pooled/group_differences_factors.m for Figure S9A; /pooled/factor_analysis_final.R for Figure S9B-C

## Install time and Run time
Install time should be short (minutes) if Matlab and RStudio are already installed - only installation step needed is to download the data and code files and save them on a local directory, keeping the folder structure the same as on the repository. Run time on a "normal" desktop computer will be long, in the order of hours to days, especially for the model-fitting and parameter recovery analyses.

