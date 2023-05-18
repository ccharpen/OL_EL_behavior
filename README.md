# OL_EL_behavior
Data and code for study investigating heterogeneity in arbitration behavior between experiential and observational learning, and its relevance for psychopathology. Preprint available here: https://psyarxiv.com/pcjg7/

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



