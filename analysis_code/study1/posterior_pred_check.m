clear all
close all
fs = filesep;

data_dir = ['..' fs '..' fs 'data'];
mod_dir = ['..' fs 'dependencies' fs 'model_functions']; %directory where modelling functions are saved (common to both studies)

addpath(['..' fs 'dependencies' fs 'cbm-master' fs 'codes']);
addpath(['..' fs 'dependencies' fs 'plotSpread']);
addpath(mod_dir); 

%load data and format data as needed (in cell array)
data = readtable([data_dir fs 'data_study1.csv']);
subID_list = unique(data.subNb);
n_all = length(subID_list);

%% Reproduce Fig 2A: Learning curves from OL, EL, DynArb model
load(['model_fitting_outputs' fs 'Recap_model_fitting.mat']);
learn_sub_OL = nan(n_all,8);
learn_sub_EL = nan(n_all,8);
learn_sub_DynArb = nan(n_all,8);

for s=1:n_all
    params_OL = fitRecap.paramRaw.ObsLearn(s,:);
    params_EL = fitRecap.paramRaw.ExpLearn(s,:);
    params_DynArb = fitRecap.paramRaw.DynArb(s,:);

    subNb = subID_list(s);
    P = table2array(data(data.subNb==subNb,2:end));
    
    %run 1000 iterations to account for stochasticity in choice-generation process
    learn_sim_OL = nan(1000,8);
    learn_sim_EL = nan(1000,8);
    
    for i=1:1000
    
        %predict choice from OL model
        P_pred_OL = generate_choice_ObsLearn(params_OL, P);
        pred_corr_OL = P_pred_OL(:,6); % whether choice is correct
        
        %predict choice from EL model
        P_pred_EL = generate_choice_ExpLearn(params_EL, P);
        pred_corr_EL = P_pred_EL(:,6); % whether choice is correct
        
        % predict choice from DynArb model
        P_pred_DynArb = generate_choice_DynArb(params_DynArb, P);
        pred_corr_DynArb = P_pred_DynArb(:,6); % whether choice is correct
       
        % Learning curve
        % identify goal token switch point - count the block break as a break
        switch_idx = [find(diff([0; (P(:,1)*2-2)+P(:,8)]));size(P,1)+1];
        learn_mat_OL = NaN(length(switch_idx)-1,8);
        learn_mat_EL = NaN(length(switch_idx)-1,8);
        learn_mat_DynArb = NaN(length(switch_idx)-1,8);

        % only calculate first 8 trials
        for j=1:(length(switch_idx)-1)
            goal_window = [switch_idx(j):switch_idx(j)+7];
            learn_mat_OL(j,:) = pred_corr_OL(goal_window);
            learn_mat_EL(j,:) = pred_corr_EL(goal_window);
            learn_mat_DynArb(j,:) = pred_corr_DynArb(goal_window);
        end
        % accuracy for each simulation
        learn_sim_OL(i,:) = mean(learn_mat_OL,1);
        learn_sim_EL(i,:) = mean(learn_mat_EL,1);
        learn_sim_DynArb(i,:) = mean(learn_mat_DynArb,1);

    end
    % average for each subject
    learn_sub_OL(s,:) = mean(learn_sim_OL,1);
    learn_sub_EL(s,:) = mean(learn_sim_EL,1);
    learn_sub_DynArb(s,:) = mean(learn_sim_DynArb,1);
    
end

% Plot learning curve
% OL
mean_learn_OL = mean(learn_sub_OL,1);
se_learn_OL = std(learn_sub_OL,1)./sqrt(size(learn_sub_OL,1));
figure;
errorbar([1:8],mean_learn_OL,se_learn_OL);
xlim([0 9]);
xlabel('# trial after reversal');
ylabel('Accuracy');
title('OL Learning curve');
% EL
mean_learn_EL = mean(learn_sub_EL,1);
se_learn_EL = std(learn_sub_EL,1)./sqrt(size(learn_sub_EL,1));
figure;
errorbar([1:8],mean_learn_EL,se_learn_EL);
xlim([0 9]);
xlabel('# trial after reversal');
ylabel('Accuracy');
title('EL Learning curve');
% DynArb
mean_learn_DynArb = mean(learn_sub_DynArb,1);
se_learn_DynArb = std(learn_sub_DynArb,1)./sqrt(size(learn_sub_DynArb,1));
figure;
errorbar([1:8],mean_learn_DynArb,se_learn_DynArb);
xlim([0 9]);
xlabel('# trial after reversal');
ylabel('Accuracy');
title('DynArb Learning curve');



%% Reproduce Fig 2B: check that OL vs EL models predict corresponding behavioral signature
ntr = 160;
nt_diff_from_OL_sub = nan(n_all,1);
prop_OL_ch_from_OL_sub = nan(n_all,1);
nt_diff_from_EL_sub = nan(n_all,1);
prop_OL_ch_from_EL_sub = nan(n_all,1);
for s=1:n_all
    params_OL = fitRecap.paramRaw.ObsLearn(s,:);
    params_EL = fitRecap.paramRaw.ExpLearn(s,:);
    subNb = subID_list(s);
    P = table2array(data(data.subNb==subNb,2:end));
    
    %run 1000 iterations to account for stochasticity in choice-generation process
    nt_diff_from_OL = nan(1000,1);
    nt_diff_from_EL = nan(1000,1);
    prop_OL_ch_from_OL = nan(1000,1);
    prop_OL_ch_from_EL = nan(1000,1);
    
    for i=1:1000
    
        %predict choice from OL model
        P_pred_OL = generate_choice_ObsLearn(params_OL, P);
        pred_ch_OL = P_pred_OL(:,4);
        
        %predict choice from EL model
        P_pred_EL = generate_choice_ExpLearn(params_EL, P);
        pred_token = P_pred_EL(:,5);
        pred_outc = P_pred_EL(:,7);
        pred_ch_EL = P_pred_EL(:,4);
        
        %define behavioral signature of OL and EL:
        %OL: whether partner's most recent choice led to an orange (1) or blue (-1) token
        past_act = nan(ntr,1);
        %EL: whether most recent token was rewarded (1) or not (-1)
        past_tok = nan(ntr,1);   
        for t=1:ntr
            %calculate whether partner's most recent choice led to orange (1)
            %or blue (-1) token and whether the most recent token was rewarded
            %or not + whether subject choice is consistent with EL and/or OL
            if P(t,4)==1
                if P(t,11)==1 %orange token obtained by partner
                    past_act(t)=1;
                elseif P(t,11)==2 %blue token obtained by partner
                    past_act(t)=-1;
                end
            else
                if (P(t-1,11)==1 && P(t,9)==P(t-1,9)) || (P(t-1,11)==2 && P(t,9)~=P(t-1,9))
                    %if partner got orange token on previous trial and repeats same choice
                    %or if partner got blue token on previous trial and switched choice
                    past_act(t)=1;
                elseif (P(t-1,11)==2 && P(t,9)==P(t-1,9)) || (P(t-1,11)==1 && P(t,9)~=P(t-1,9))
                    %if partner got blue token on previous trial and repeats same choice
                    %or if partner got orange token on previous trial and switched choice
                    past_act(t)=-1;
                end
            end
            if t==1 %no past token evidence on trial 1
                past_tok(t) = 0; 
            else
                if (pred_token(t-1)==1 && pred_outc(t-1)>0) || (pred_token(t-1)==2 && pred_outc(t-1)==0)
                    %past orange token was rewarded or past blue token was not
                    past_tok(t)=1;
                elseif (pred_token(t-1)==2 && pred_outc(t-1)>0) || (pred_token(t-1)==1 && pred_outc(t-1)==0)
                    %past blue token was rewarded or past orange token was not
                    past_tok(t)=-1;
                end
            end
        end
    
        %define trials consistent with OL and with EL from OL-predicted choices   
        ch_OL_from_OL = (past_act==1 & past_tok==-1 & pred_ch_OL==1) | (past_act==-1 & past_tok==1 & pred_ch_OL==0);
        ch_EL_from_OL = (past_act==1 & past_tok==-1 & pred_ch_OL==0) | (past_act==-1 & past_tok==1 & pred_ch_OL==1);
        ind_diff_from_OL = ch_OL_from_OL | ch_EL_from_OL; %trials where OL and EL make different predictions
        nt_diff_from_OL(i,1) = sum(ind_diff_from_OL);
        prop_OL_ch_from_OL(i,1) = mean(ch_OL_from_OL(ind_diff_from_OL));

        %define trials consistent with OL and with EL from EL-predicted choices   
        ch_OL_from_EL = (past_act==1 & past_tok==-1 & pred_ch_EL==1) | (past_act==-1 & past_tok==1 & pred_ch_EL==0);
        ch_EL_from_EL = (past_act==1 & past_tok==-1 & pred_ch_EL==0) | (past_act==-1 & past_tok==1 & pred_ch_EL==1);
        ind_diff_from_EL = ch_OL_from_EL | ch_EL_from_EL; %trials where OL and EL make different predictions
        nt_diff_from_EL(i,1) = sum(ind_diff_from_EL);
        prop_OL_ch_from_EL(i,1) = mean(ch_OL_from_EL(ind_diff_from_EL));
    end
    nt_diff_from_OL_sub(s,1) = mean(nt_diff_from_OL);
    prop_OL_ch_from_OL_sub(s,1) = mean(prop_OL_ch_from_OL);
    nt_diff_from_EL_sub(s,1) = mean(nt_diff_from_EL);
    prop_OL_ch_from_EL_sub(s,1) = mean(prop_OL_ch_from_EL);
end
%load Behavioral variables to compare with model predictions
load('Behavioral_variables.mat')

%plot histograms of model predictions and data
figure;
subplot(2,1,1); hold on
histogram(prop_OL_ch_from_EL_sub,20,'FaceAlpha',0.6);
histogram(prop_OL_ch_from_OL_sub,15,'FaceAlpha',0.6);
set(gca,'box','off')
xlim([0.2 0.8])
ylabel('count')
l = legend({'EL model','OL model'});
title(l,'predicted by:')
title('Model predictions')
subplot(2,1,2);
histogram(Behavior.Prop_OL_ch(:,1),25,'FaceColor',[0.5 0.5 0.5]);
set(gca,'box','off')
xlim([0.2 0.8])
ylabel('count')
xlabel('proportion of choices consistent with OL (vs EL)')
title('Data')

%plot correlations between data and model predictions
iEL = Behavior.Prop_OL_ch(:,1)<0.5;
iOL = Behavior.Prop_OL_ch(:,1)>0.5;
figure;
subplot(1,2,1); hold on
plot(Behavior.Prop_OL_ch(iEL,1),prop_OL_ch_from_EL_sub(iEL),'.','Color','#0072BD'); lsline()
plot(Behavior.Prop_OL_ch(iEL,1),prop_OL_ch_from_OL_sub(iEL),'.','Color','#D95319'); lsline()
xlim([0.2 0.5]); ylim([0.2 0.8])
title({'Experiential learners'; '(OL choice prop. < 0.5)'})
xlabel('Data')
ylabel('Model predictions')
legend({'EL model','','OL model',''})
subplot(1,2,2); hold on
plot(Behavior.Prop_OL_ch(iOL,1),prop_OL_ch_from_EL_sub(iOL),'.','Color','#0072BD'); lsline()
plot(Behavior.Prop_OL_ch(iOL,1),prop_OL_ch_from_OL_sub(iOL),'.','Color','#D95319'); lsline()
xlim([0.5 0.8]); ylim([0.2 0.8])
title({'Observational learners'; '(OL choice prop. > 0.5)'})
xlabel('Data')

%print correlations
[r,p] = corr(Behavior.Prop_OL_ch(iEL,1),prop_OL_ch_from_EL_sub(iEL))
[r,p] = corr(Behavior.Prop_OL_ch(iEL,1),prop_OL_ch_from_OL_sub(iEL))
[r,p] = corr(Behavior.Prop_OL_ch(iOL,1),prop_OL_ch_from_EL_sub(iOL))
[r,p] = corr(Behavior.Prop_OL_ch(iOL,1),prop_OL_ch_from_OL_sub(iOL))

%% Run mixed effect glms - trial definition of uncertainties
clear all
close all

sim_cond = 'DynArb'; % which model for simulation: OL, EL, DynArb

fs = filesep;
load(['model_fitting_outputs' fs 'Recap_model_fitting.mat']);
dir_data = ['..' fs '..' fs 'data'];
addpath(['..' fs 'dependencies' fs 'plotSpread'])

% init the randomization screen
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

%load data
data = readtable([dir_data fs 'data_study1.csv']);

nsub = length(unique(data.subNb)); %number of subjects
ntr = 160; %number of trials
behavior_all = {};


Behavior = struct();
Behavior.subID_list = unique(data.subNb);
glm_mat_allsubs = [];
glm_mat_OLunc = [];
glm_mat_ELunc = [];
glm_mat_Mag = [];

for sub=1:nsub

    %load subject data
    subNb = Behavior.subID_list(sub);
    P = table2array(data(data.subNb==subNb,2:end));

    % OL prediction
    if strcmp(sim_cond,'OL')
        params_OL = fitRecap.paramRaw.ObsLearn(sub,:);
        P_pred = generate_choice_ObsLearn(params_OL, P);
    elseif strcmp(sim_cond,'EL')
        params_EL = fitRecap.paramRaw.ExpLearn(sub,:);
        P_pred = generate_choice_ExpLearn(params_EL, P);
    elseif strcmp(sim_cond,'DynArb')
        params_DynArb = fitRecap.paramRaw.DynArb(sub,:);
        P_pred = generate_choice_DynArb(params_DynArb, P);
    end

    %define uncertainty condition based on key trials, instead of design conditions
    %OL unc = low if past 2 box-token transition consistent; high if not
    %EL unc = low if W-st-W, W-sw-0, 0-st-0, 0-sw-W; high otherwise
    %RM = high if >25, low otherwise
    OL_lu = zeros(ntr,1);
    EL_lu = zeros(ntr,1);
    RM_h = zeros(ntr,1);
    outc_bin = P_pred(:,7);
    for t=1:ntr
        if P(t,4)>=3
            if P(t-1,6)==P(t,6)
                OL_lu(t) = 1;
            end
            if (outc_bin(t-2)==outc_bin(t-1) && P_pred(t-2,4)==P_pred(t-1,4)) || (outc_bin(t-2)~=outc_bin(t-1) && P_pred(t-2,4)~=P_pred(t-1,4))
                %same choice, same outcome OR different choice, different outcome
                EL_lu(t) = 1;
            end
        end
        if P(t,3)>=2 && P_pred(t-1,8)>25/100
            RM_h(t) = 1;
        end
    end
    OL_lu = logical(OL_lu);
    EL_lu = logical(EL_lu);
    RM_h = logical(RM_h);


    % Calculate variables for the GLM and for behavioral signatures of the two strategies
    %1) whether partner's most recent choice led to an orange (1) or blue (-1) token
    past_act = nan(ntr,1);
    %2) whether most recent token was rewarded (1) or not (-1)
    past_tok = nan(ntr,1);
    %3) whether most recent token was rewarded (1) or not (-1) * outcome magnitude
    past_tok_comb = nan(ntr,1);

    for t=1:ntr
        %calculate whether partner's most recent choice led to orange (1)
        %or blue (-1) token and whether the most recent token was rewarded
        %or not + whether subject choice is consistent with EL and/or OL
        if P(t,4)==1
            if P(t,11)==1 %orange token obtained by partner
                past_act(t)=1;
            elseif P(t,11)==2 %blue token obtained by partner
                past_act(t)=-1;
            end
        else
            if (P(t-1,11)==1 && P(t,9)==P(t-1,9)) || (P(t-1,11)==2 && P(t,9)~=P(t-1,9))
                %if partner got orange token on previous trial and repeats same choice
                %or if partner got blue token on previous trial and switched choice
                past_act(t)=1;
            elseif P(t-1,11)==2 && P(t,9)==P(t-1,9) || (P(t-1,11)==1 && P(t,9)~=P(t-1,9))
                %if partner got blue token on previous trial and repeats same choice
                %or if partner got orange token on previous trial and switched choice
                past_act(t)=-1;
            end
        end

        if t==1 %no past token evidence on trial 1
            past_tok(t) = 0; 
            past_tok_comb(t) = 0;
        else
            if (P_pred(t-1,5)==1 && P_pred(t-1,8)>0) || (P_pred(t-1,5)==2 && P_pred(t-1,8)==0)
                %past orange token was rewarded or past blue token was not
                past_tok(t)=1;
            elseif (P_pred(t-1,5)==2 && P_pred(t-1,8)>0) || (P_pred(t-1,5)==1 && P_pred(t-1,8)==0)
                %past blue token was rewarded or past orange token was not
                past_tok(t)=-1;
            end
            if P_pred(t-1,5)==1 %past token was orange
                past_tok_comb(t) = P_pred(t-1,8); % outcome magnitude, not binary outcome
            elseif P_pred(t-1,5)==2 %past token was blue
                past_tok_comb(t) = -P_pred(t-1,8);
            end
        end
    end

    % build matrices for glms
    good = ~isnan(past_tok) & ~isnan(P_pred(:,4)); 

    %glm matrix for assessing hybrid behavior
    glm_mat1 = [zscore(past_act(good)) zscore(past_tok_comb(good)) P_pred(good,4)];
    glm_mat_allsubs = [glm_mat_allsubs; [ones(sum(good),1)*subNb glm_mat1]];

    %build glm matrix to look at interactions with uncertainty in mixed-effect GLM
    ELUncL = double(EL_lu); ELUncL(ELUncL == 0) = -1;
    OLUncL = double(OL_lu); OLUncL(OLUncL == 0) = -1;
    MagH = double(RM_h); MagH(MagH == 0) = -1;

    glm_mat_OLunc = [glm_mat_OLunc; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & OL_lu)); zeros(sum(good & ~OL_lu),1)] ...
        [zeros(sum(good & OL_lu),1); zscore(past_act(good & ~OL_lu))] ...
        [zscore(past_tok_comb(good & OL_lu)); zeros(sum(good & ~OL_lu),1)] ...
        [zeros(sum(good & OL_lu),1); zscore(past_tok_comb(good & ~OL_lu))] ...
        [P_pred(good & OL_lu,4);P_pred(good & ~OL_lu,4)]];

    glm_mat_ELunc = [glm_mat_ELunc; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
        [zeros(sum(good & EL_lu),1); zscore(past_act(good & ~EL_lu))] ...
        [zscore(past_tok_comb(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
        [zeros(sum(good & EL_lu),1); zscore(past_tok_comb(good & ~EL_lu))] ...
        [P_pred(good & EL_lu,4);P_pred(good & ~EL_lu,4)]];

    glm_mat_Mag = [glm_mat_Mag; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & RM_h)); zeros(sum(good & ~RM_h),1)] ...
        [zeros(sum(good & RM_h),1); zscore(past_act(good & ~RM_h))] ...
        [zscore(past_tok_comb(good & RM_h)); zeros(sum(good & ~RM_h),1)] ...
        [zeros(sum(good & RM_h),1); zscore(past_tok_comb(good & ~RM_h))] ...
        [P_pred(good & RM_h,4);P_pred(good & ~RM_h,4)]];

end

%hybrid behavior
glm_tbl = array2table(glm_mat_allsubs, 'VariableNames',{'subNb','PastAct','PastTok','choice'});
glme = fitglme(glm_tbl,'choice ~ 1 + PastAct + PastTok + (1 + PastAct + PastTok|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
% disp(glme)
% anova(glme)
[~,~,stats] = randomEffects(glme);
Meffects = dataset2table(glme.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct = U.PastAct + Meffects.Estimate(2);
U.PastTok = U.PastTok + Meffects.Estimate(3);
Behavior.GLME.input_data = glm_tbl;
Behavior.GLME.Coefficients = Meffects;
Behavior.GLME.Reffects = U;
Behavior.GLME.R_stats = Reffects;

%separate GLMs for each uncertainty condition
%OL uncertainty
glm_OLunc_tbl = array2table(glm_mat_OLunc, 'VariableNames',...
    {'subNb','PastAct_LowOLU','PastAct_HighOLU','PastTok_LowOLU','PastTok_HighOLU','choice'});
glme_OLunc = fitglme(glm_OLunc_tbl,...
    'choice ~ 1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU + (1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
%disp(glme_OLunc)
[~,~,stats] = randomEffects(glme_OLunc);
Meffects = dataset2table(glme_OLunc.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct_LowOLU = U.PastAct_LowOLU + Meffects.Estimate(2);
U.PastAct_HighOLU = U.PastAct_HighOLU + Meffects.Estimate(3);
U.PastTok_LowOLU = U.PastTok_LowOLU + Meffects.Estimate(4);
U.PastTok_HighOLU = U.PastTok_HighOLU + Meffects.Estimate(5);
Behavior.GLME_OLunc.input_data = glm_OLunc_tbl;
Behavior.GLME_OLunc.Coefficients = Meffects;
Behavior.GLME_OLunc.Reffects = U;
Behavior.GLME_OLunc.R_stats = Reffects;

%EL uncertainty
glm_ELunc_tbl = array2table(glm_mat_ELunc, 'VariableNames',...
    {'subNb','PastAct_LowELU','PastAct_HighELU','PastTok_LowELU','PastTok_HighELU','choice'});
glme_ELunc = fitglme(glm_ELunc_tbl,...
    'choice ~ 1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU + (1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
%disp(glme_ELunc)
[~,~,stats] = randomEffects(glme_ELunc);
Meffects = dataset2table(glme_ELunc.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
%U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct_LowELU = U.PastAct_LowELU + Meffects.Estimate(2);
U.PastAct_HighELU = U.PastAct_HighELU + Meffects.Estimate(3);
U.PastTok_LowELU = U.PastTok_LowELU + Meffects.Estimate(4);
U.PastTok_HighELU = U.PastTok_HighELU + Meffects.Estimate(5);
Behavior.GLME_ELunc.input_data = glm_ELunc_tbl;
Behavior.GLME_ELunc.Coefficients = Meffects;
Behavior.GLME_ELunc.Reffects = U;
Behavior.GLME_ELunc.R_stats = Reffects;

%reward magnitude
glm_Mag_tbl = array2table(glm_mat_Mag, 'VariableNames',...
    {'subNb','PastAct_HighMag','PastAct_LowMag','PastTok_HighMag','PastTok_LowMag','choice'});
glme_Mag = fitglme(glm_Mag_tbl,...
    'choice ~ 1 + PastAct_HighMag + PastAct_LowMag + PastTok_HighMag + PastTok_LowMag + (1 + PastAct_HighMag + PastAct_LowMag + PastTok_HighMag + PastTok_LowMag|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
%disp(glme_Mag)
[~,~,stats] = randomEffects(glme_Mag);
Meffects = dataset2table(glme_Mag.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct_HighMag = U.PastAct_HighMag + Meffects.Estimate(2);
U.PastAct_LowMag = U.PastAct_LowMag + Meffects.Estimate(3);
U.PastTok_HighMag = U.PastTok_HighMag + Meffects.Estimate(4);
U.PastTok_LowMag = U.PastTok_LowMag + Meffects.Estimate(5);
Behavior.GLME_Mag.input_data = glm_Mag_tbl;
Behavior.GLME_Mag.Coefficients = Meffects;
Behavior.GLME_Mag.Reffects = U;
Behavior.GLME_Mag.R_stats = Reffects;


% plot
%GLM
figure; 
subplot(1,4,1)
hold on
bar(1,Behavior.GLME.Coefficients.Estimate([3 2]),0.57,'EdgeColor','k','LineWidth',1); 
plotSpread(Behavior.GLME.Reffects.PastTok,'xValues',0.85,'distributionColors',[0 0.32 0.47],'spreadWidth',0.25);
plotSpread(Behavior.GLME.Reffects.PastAct,'xValues',1.15,'distributionColors',[0.52 0.26 0.08],'spreadWidth',0.25);
errorbar([0.85 1.15],Behavior.GLME.Coefficients.Estimate([3 2]),Behavior.GLME.Coefficients.SE([3 2]),'.k','LineWidth',1.5);
xticks([])
ylabel('ME-GLM effect')
legend({'Effect of past outcome (EL)','Effect of past partner''s action (OL)'})

%GLM OL and EL uncertainty
subplot(1,4,2)
hold on
h(1) = plot(0.85:1.85,Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),'-', 'LineWidth', 1.5);
h(2) = plot(1.15:2.15,Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),'-', 'LineWidth', 1.5);
plotSpread([Behavior.GLME_OLunc.Reffects.PastTok_LowOLU Behavior.GLME_OLunc.Reffects.PastTok_HighOLU],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_OLunc.Reffects.PastAct_LowOLU Behavior.GLME_OLunc.Reffects.PastAct_HighOLU],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),Behavior.GLME_OLunc.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),Behavior.GLME_OLunc.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35])
xticks([1 2])
xticklabels({'Low','High'})
xlabel('OL uncertainty')
ylabel('ME-GLM effect')
subplot(1,4,3)
hold on
h(1) = plot(0.85:1.85,Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5);
h(2) = plot(1.15:2.15,Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_ELunc.Reffects.PastTok_LowELU Behavior.GLME_ELunc.Reffects.PastTok_HighELU],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_ELunc.Reffects.PastAct_LowELU Behavior.GLME_ELunc.Reffects.PastAct_HighELU],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),Behavior.GLME_ELunc.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),Behavior.GLME_ELunc.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35])
xticks([1 2])
xticklabels({'Low','High'})
xlabel('EL uncertainty')
ylabel('ME-GLM effect')
%legend(h, {'Past outcome (EL effect)','Past action (OL effect)'});

%GLM reward magnitude 
subplot(1,4,4)
h(1) = plot(0.85:1.85,Behavior.GLME_Mag.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5); hold
h(2) = plot(1.15:2.15,Behavior.GLME_Mag.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_Mag.Reffects.PastTok_HighMag Behavior.GLME_Mag.Reffects.PastTok_LowMag],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_Mag.Reffects.PastAct_HighMag Behavior.GLME_Mag.Reffects.PastAct_LowMag],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_Mag.Coefficients.Estimate([4 5]),Behavior.GLME_Mag.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_Mag.Coefficients.Estimate([2 3]),Behavior.GLME_Mag.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35])
xticks([1 2])
xticklabels({'High','Low'})
xlabel('Reward Magnitude')
ylabel('ME-GLM effect')
legend(h, {'Past outcome (EL effect)','Past action (OL effect)'});


% correlation with actual data
Behavior_actual = load('Behavioral_variables.mat');
% general GLME effect
figure;
subplot(1,5,1)
hold on
scatter(Behavior_actual.Behavior.GLME.Reffects.PastAct,Behavior.GLME.Reffects.PastAct,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('PastAct');
refline([1,0]);
subplot(1,5,2)
hold on
scatter(Behavior_actual.Behavior.GLME.Reffects.PastTok,Behavior.GLME.Reffects.PastTok,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('PastTok');
refline([1,0]);
% diff of PastAct between low and high OL unc
diff_OL_actual = Behavior_actual.Behavior.GLME_OLunc.Reffects.PastAct_LowOLU-Behavior_actual.Behavior.GLME_OLunc.Reffects.PastAct_HighOLU;
diff_OL_sim = Behavior.GLME_OLunc.Reffects.PastAct_LowOLU-Behavior.GLME_OLunc.Reffects.PastAct_HighOLU;
subplot(1,5,3)
hold on
scatter(diff_OL_actual,diff_OL_sim,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('OL effect on PastAct');
refline([1,0]);
% diff of PastTok between low and high EL unc
diff_EL_actual = Behavior_actual.Behavior.GLME_ELunc.Reffects.PastTok_LowELU-Behavior_actual.Behavior.GLME_ELunc.Reffects.PastTok_HighELU;
diff_EL_sim = Behavior.GLME_ELunc.Reffects.PastTok_LowELU-Behavior.GLME_ELunc.Reffects.PastTok_HighELU;
subplot(1,5,4)
hold on
scatter(diff_EL_actual,diff_EL_sim,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('EL effect on PastTok');
refline([1,0]);
% diff of PastTok between low and high Rew mag
diff_rew_actual = Behavior_actual.Behavior.GLME_Mag.Reffects.PastTok_LowMag-Behavior_actual.Behavior.GLME_Mag.Reffects.PastTok_HighMag;
diff_rew_sim = Behavior.GLME_Mag.Reffects.PastTok_LowMag-Behavior.GLME_Mag.Reffects.PastTok_HighMag;
subplot(1,5,5)
hold on
scatter(diff_rew_actual,diff_rew_sim,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('Magnitude effect on PastTok');
refline([1,0]);

%% Estimate actual mixed effect glms - block definition of uncertainties
clear all
close all

fs = filesep;
load(['model_fitting_outputs' fs 'Recap_model_fitting.mat']);
dir_data = ['..' fs '..' fs 'data'];
addpath(['..' fs 'dependencies' fs 'plotSpread'])

% init the randomization screen
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

%load data
data = readtable([dir_data fs 'data_study1.csv']);

nsub = length(unique(data.subNb)); %number of subjects
ntr = 160; %number of trials
behavior_all = {};


Behavior = struct();
Behavior.subID_list = unique(data.subNb);
glm_mat_allsubs = [];
glm_mat_OLunc = [];
glm_mat_ELunc = [];

for sub=1:nsub

    %load subject data
    subNb = Behavior.subID_list(sub);
    P = table2array(data(data.subNb==subNb,2:end));

    %define uncertainty condition based on design conditions, instead of key trials
    %OL unc = P(:,5) 0.6 = hu, 0.8 = lu
    %EL unc = P(:,7) 0.6/0.4 = hu, 0.2/0.8 = lu
    %RM = high if >25, low otherwise
    OL_lu = (P(:,5) == 0.8);
    EL_lu = (P(:,7) == 0.2 | P(:,7) == 0.8);
    outc_bin = P(:,18); outc_bin(outc_bin>0)=1;

    % Calculate variables for the GLM and for behavioral signatures of the two strategies
    %1) whether partner's most recent choice led to an orange (1) or blue (-1) token
    past_act = nan(ntr,1);
    %2) whether most recent token was rewarded (1) or not (-1)
    past_tok = nan(ntr,1);
    %3) whether most recent token was rewarded (1) or not (-1) * outcome magnitude
    past_tok_comb = nan(ntr,1);
    
    for t=1:ntr
        %calculate whether partner's most recent choice led to orange (1)
        %or blue (-1) token and whether the most recent token was rewarded
        %or not + whether subject choice is consistent with EL and/or OL
        if P(t,4)==1
            if P(t,11)==1 %orange token obtained by partner
                past_act(t)=1;
            elseif P(t,11)==2 %blue token obtained by partner
                past_act(t)=-1;
            end
        else
            if (P(t-1,11)==1 && P(t,9)==P(t-1,9)) || (P(t-1,11)==2 && P(t,9)~=P(t-1,9))
                %if partner got orange token on previous trial and repeats same choice
                %or if partner got blue token on previous trial and switched choice
                past_act(t)=1;
            elseif P(t-1,11)==2 && P(t,9)==P(t-1,9) || (P(t-1,11)==1 && P(t,9)~=P(t-1,9))
                %if partner got blue token on previous trial and repeats same choice
                %or if partner got orange token on previous trial and switched choice
                past_act(t)=-1;
            end
        end
        
        if t==1 %no past token evidence on trial 1
            past_tok(t) = 0; 
            past_tok_comb(t) = 0;
        else
            if (P(t-1,17)==1 && P(t-1,18)>0) || (P(t-1,17)==2 && P(t-1,18)==0)
                %past orange token was rewarded or past blue token was not
                past_tok(t)=1;
            elseif (P(t-1,17)==2 && P(t-1,18)>0) || (P(t-1,17)==1 && P(t-1,18)==0)
                %past blue token was rewarded or past orange token was not
                past_tok(t)=-1;
            end
            if P(t-1,17)==1 %past token was orange
                past_tok_comb(t) = P(t-1,18);
            elseif P(t-1,17)==2 %past token was blue
                past_tok_comb(t) = -P(t-1,18);
            end
        end
    end

    % build matrices for glms
    good = ~isnan(past_tok) & ~isnan(P(:,14)); 

    %glm matrix for assessing hybrid behavior
    glm_mat1 = [zscore(past_act(good)) zscore(past_tok_comb(good)) P(good,14)];
    glm_mat_allsubs = [glm_mat_allsubs; [ones(sum(good),1)*subNb glm_mat1]];

    %build glm matrix to look at interactions with uncertainty in mixed-effect GLM
    ELUncL = double(EL_lu); ELUncL(ELUncL == 0) = -1;
    OLUncL = double(OL_lu); OLUncL(OLUncL == 0) = -1;

    glm_mat_OLunc = [glm_mat_OLunc; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & OL_lu)); zeros(sum(good & ~OL_lu),1)] ...
        [zeros(sum(good & OL_lu),1); zscore(past_act(good & ~OL_lu))] ...
        [zscore(past_tok_comb(good & OL_lu)); zeros(sum(good & ~OL_lu),1)] ...
        [zeros(sum(good & OL_lu),1); zscore(past_tok_comb(good & ~OL_lu))] ...
        [P(good & OL_lu,14);P(good & ~OL_lu,14)]];
    
    glm_mat_ELunc = [glm_mat_ELunc; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
        [zeros(sum(good & EL_lu),1); zscore(past_act(good & ~EL_lu))] ...
        [zscore(past_tok_comb(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
        [zeros(sum(good & EL_lu),1); zscore(past_tok_comb(good & ~EL_lu))] ...
        [P(good & EL_lu,14);P(good & ~EL_lu,14)]];
end

%hybrid behavior
glm_tbl = array2table(glm_mat_allsubs, 'VariableNames',{'subNb','PastAct','PastTok','choice'});
glme = fitglme(glm_tbl,'choice ~ 1 + PastAct + PastTok + (1 + PastAct + PastTok|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
% disp(glme)
% anova(glme)
[~,~,stats] = randomEffects(glme);
Meffects = dataset2table(glme.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct = U.PastAct + Meffects.Estimate(2);
U.PastTok = U.PastTok + Meffects.Estimate(3);
Behavior.GLME.input_data = glm_tbl;
Behavior.GLME.Coefficients = Meffects;
Behavior.GLME.Reffects = U;
Behavior.GLME.R_stats = Reffects;

%separate GLMs for each uncertainty condition
%OL uncertainty
glm_OLunc_tbl = array2table(glm_mat_OLunc, 'VariableNames',...
    {'subNb','PastAct_LowOLU','PastAct_HighOLU','PastTok_LowOLU','PastTok_HighOLU','choice'});
glme_OLunc = fitglme(glm_OLunc_tbl,...
    'choice ~ 1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU + (1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
%disp(glme_OLunc)
[~,~,stats] = randomEffects(glme_OLunc);
Meffects = dataset2table(glme_OLunc.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct_LowOLU = U.PastAct_LowOLU + Meffects.Estimate(2);
U.PastAct_HighOLU = U.PastAct_HighOLU + Meffects.Estimate(3);
U.PastTok_LowOLU = U.PastTok_LowOLU + Meffects.Estimate(4);
U.PastTok_HighOLU = U.PastTok_HighOLU + Meffects.Estimate(5);
Behavior.GLME_OLunc.input_data = glm_OLunc_tbl;
Behavior.GLME_OLunc.Coefficients = Meffects;
Behavior.GLME_OLunc.Reffects = U;
Behavior.GLME_OLunc.R_stats = Reffects;

%EL uncertainty
glm_ELunc_tbl = array2table(glm_mat_ELunc, 'VariableNames',...
    {'subNb','PastAct_LowELU','PastAct_HighELU','PastTok_LowELU','PastTok_HighELU','choice'});
glme_ELunc = fitglme(glm_ELunc_tbl,...
    'choice ~ 1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU + (1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
%disp(glme_ELunc)
[~,~,stats] = randomEffects(glme_ELunc);
Meffects = dataset2table(glme_ELunc.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
%U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct_LowELU = U.PastAct_LowELU + Meffects.Estimate(2);
U.PastAct_HighELU = U.PastAct_HighELU + Meffects.Estimate(3);
U.PastTok_LowELU = U.PastTok_LowELU + Meffects.Estimate(4);
U.PastTok_HighELU = U.PastTok_HighELU + Meffects.Estimate(5);
Behavior.GLME_ELunc.input_data = glm_ELunc_tbl;
Behavior.GLME_ELunc.Coefficients = Meffects;
Behavior.GLME_ELunc.Reffects = U;
Behavior.GLME_ELunc.R_stats = Reffects;

save('Behavioral_blockunc.mat','Behavior');
%% Run simulation on mixed effect glms - block definition of uncertainties
clear all
close all

sim_cond = 'DynArb'; % which model for simulation: OL, EL, DynArb

fs = filesep;
load(['model_fitting_outputs' fs 'Recap_model_fitting.mat']);
dir_data = ['..' fs '..' fs 'data'];
addpath(['..' fs 'dependencies' fs 'plotSpread'])

% init the randomization screen
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

%load data
data = readtable([dir_data fs 'data_study1.csv']);

nsub = length(unique(data.subNb)); %number of subjects
ntr = 160; %number of trials
behavior_all = {};


Behavior = struct();
Behavior.subID_list = unique(data.subNb);
glm_mat_allsubs = [];
glm_mat_OLunc = [];
glm_mat_ELunc = [];
glm_mat_Mag = [];

for sub=1:nsub

    %load subject data
    subNb = Behavior.subID_list(sub);
    P = table2array(data(data.subNb==subNb,2:end));

    % OL prediction
    if strcmp(sim_cond,'OL')
        params_OL = fitRecap.paramRaw.ObsLearn(sub,:);
        P_pred = generate_choice_ObsLearn(params_OL, P);
    elseif strcmp(sim_cond,'EL')
        params_EL = fitRecap.paramRaw.ExpLearn(sub,:);
        P_pred = generate_choice_ExpLearn(params_EL, P);
    elseif strcmp(sim_cond,'DynArb')
        params_DynArb = fitRecap.paramRaw.DynArb(sub,:);
        P_pred = generate_choice_DynArb(params_DynArb, P);
    end

    %define uncertainty condition based on design conditions, instead of key trials
    %OL unc = P(:,5) 0.6 = hu, 0.8 = lu
    %EL unc = P(:,7) 0.6/0.4 = hu, 0.2/0.8 = lu
    %RM = high if >25, low otherwise
    OL_lu = (P(:,5) == 0.8);
    EL_lu = (P(:,7) == 0.2 | P(:,7) == 0.8);
 %   RM_h = zeros(ntr,1);
    outc_bin = P_pred(:,7);

    % Calculate variables for the GLM and for behavioral signatures of the two strategies
    %1) whether partner's most recent choice led to an orange (1) or blue (-1) token
    past_act = nan(ntr,1);
    %2) whether most recent token was rewarded (1) or not (-1)
    past_tok = nan(ntr,1);
    %3) whether most recent token was rewarded (1) or not (-1) * outcome magnitude
    past_tok_comb = nan(ntr,1);

    for t=1:ntr
        %calculate whether partner's most recent choice led to orange (1)
        %or blue (-1) token and whether the most recent token was rewarded
        %or not + whether subject choice is consistent with EL and/or OL
        if P(t,4)==1
            if P(t,11)==1 %orange token obtained by partner
                past_act(t)=1;
            elseif P(t,11)==2 %blue token obtained by partner
                past_act(t)=-1;
            end
        else
            if (P(t-1,11)==1 && P(t,9)==P(t-1,9)) || (P(t-1,11)==2 && P(t,9)~=P(t-1,9))
                %if partner got orange token on previous trial and repeats same choice
                %or if partner got blue token on previous trial and switched choice
                past_act(t)=1;
            elseif P(t-1,11)==2 && P(t,9)==P(t-1,9) || (P(t-1,11)==1 && P(t,9)~=P(t-1,9))
                %if partner got blue token on previous trial and repeats same choice
                %or if partner got orange token on previous trial and switched choice
                past_act(t)=-1;
            end
        end

        if t==1 %no past token evidence on trial 1
            past_tok(t) = 0; 
            past_tok_comb(t) = 0;
        else
            if (P_pred(t-1,5)==1 && P_pred(t-1,8)>0) || (P_pred(t-1,5)==2 && P_pred(t-1,8)==0)
                %past orange token was rewarded or past blue token was not
                past_tok(t)=1;
            elseif (P_pred(t-1,5)==2 && P_pred(t-1,8)>0) || (P_pred(t-1,5)==1 && P_pred(t-1,8)==0)
                %past blue token was rewarded or past orange token was not
                past_tok(t)=-1;
            end
            if P_pred(t-1,5)==1 %past token was orange
                past_tok_comb(t) = P_pred(t-1,8); % outcome magnitude, not binary outcome
            elseif P_pred(t-1,5)==2 %past token was blue
                past_tok_comb(t) = -P_pred(t-1,8);
            end
        end
    end

    % build matrices for glms
    good = ~isnan(past_tok) & ~isnan(P_pred(:,4)); 

    %glm matrix for assessing hybrid behavior
    glm_mat1 = [zscore(past_act(good)) zscore(past_tok_comb(good)) P_pred(good,4)];
    glm_mat_allsubs = [glm_mat_allsubs; [ones(sum(good),1)*subNb glm_mat1]];

    %build glm matrix to look at interactions with uncertainty in mixed-effect GLM
    ELUncL = double(EL_lu); ELUncL(ELUncL == 0) = -1;
    OLUncL = double(OL_lu); OLUncL(OLUncL == 0) = -1;

    glm_mat_OLunc = [glm_mat_OLunc; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & OL_lu)); zeros(sum(good & ~OL_lu),1)] ...
        [zeros(sum(good & OL_lu),1); zscore(past_act(good & ~OL_lu))] ...
        [zscore(past_tok_comb(good & OL_lu)); zeros(sum(good & ~OL_lu),1)] ...
        [zeros(sum(good & OL_lu),1); zscore(past_tok_comb(good & ~OL_lu))] ...
        [P_pred(good & OL_lu,4);P_pred(good & ~OL_lu,4)]];

    glm_mat_ELunc = [glm_mat_ELunc; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
        [zeros(sum(good & EL_lu),1); zscore(past_act(good & ~EL_lu))] ...
        [zscore(past_tok_comb(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
        [zeros(sum(good & EL_lu),1); zscore(past_tok_comb(good & ~EL_lu))] ...
        [P_pred(good & EL_lu,4);P_pred(good & ~EL_lu,4)]];

end

%hybrid behavior
glm_tbl = array2table(glm_mat_allsubs, 'VariableNames',{'subNb','PastAct','PastTok','choice'});
glme = fitglme(glm_tbl,'choice ~ 1 + PastAct + PastTok + (1 + PastAct + PastTok|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
% disp(glme)
% anova(glme)
[~,~,stats] = randomEffects(glme);
Meffects = dataset2table(glme.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct = U.PastAct + Meffects.Estimate(2);
U.PastTok = U.PastTok + Meffects.Estimate(3);
Behavior.GLME.input_data = glm_tbl;
Behavior.GLME.Coefficients = Meffects;
Behavior.GLME.Reffects = U;
Behavior.GLME.R_stats = Reffects;

%separate GLMs for each uncertainty condition
%OL uncertainty
glm_OLunc_tbl = array2table(glm_mat_OLunc, 'VariableNames',...
    {'subNb','PastAct_LowOLU','PastAct_HighOLU','PastTok_LowOLU','PastTok_HighOLU','choice'});
glme_OLunc = fitglme(glm_OLunc_tbl,...
    'choice ~ 1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU + (1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
%disp(glme_OLunc)
[~,~,stats] = randomEffects(glme_OLunc);
Meffects = dataset2table(glme_OLunc.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct_LowOLU = U.PastAct_LowOLU + Meffects.Estimate(2);
U.PastAct_HighOLU = U.PastAct_HighOLU + Meffects.Estimate(3);
U.PastTok_LowOLU = U.PastTok_LowOLU + Meffects.Estimate(4);
U.PastTok_HighOLU = U.PastTok_HighOLU + Meffects.Estimate(5);
Behavior.GLME_OLunc.input_data = glm_OLunc_tbl;
Behavior.GLME_OLunc.Coefficients = Meffects;
Behavior.GLME_OLunc.Reffects = U;
Behavior.GLME_OLunc.R_stats = Reffects;

%EL uncertainty
glm_ELunc_tbl = array2table(glm_mat_ELunc, 'VariableNames',...
    {'subNb','PastAct_LowELU','PastAct_HighELU','PastTok_LowELU','PastTok_HighELU','choice'});
glme_ELunc = fitglme(glm_ELunc_tbl,...
    'choice ~ 1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU + (1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
%disp(glme_ELunc)
[~,~,stats] = randomEffects(glme_ELunc);
Meffects = dataset2table(glme_ELunc.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
%U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct_LowELU = U.PastAct_LowELU + Meffects.Estimate(2);
U.PastAct_HighELU = U.PastAct_HighELU + Meffects.Estimate(3);
U.PastTok_LowELU = U.PastTok_LowELU + Meffects.Estimate(4);
U.PastTok_HighELU = U.PastTok_HighELU + Meffects.Estimate(5);
Behavior.GLME_ELunc.input_data = glm_ELunc_tbl;
Behavior.GLME_ELunc.Coefficients = Meffects;
Behavior.GLME_ELunc.Reffects = U;
Behavior.GLME_ELunc.R_stats = Reffects;


% plot
%GLM
figure; 
subplot(1,3,1)
hold on
bar(1,Behavior.GLME.Coefficients.Estimate([3 2]),0.57,'EdgeColor','k','LineWidth',1); 
plotSpread(Behavior.GLME.Reffects.PastTok,'xValues',0.85,'distributionColors',[0 0.32 0.47],'spreadWidth',0.25);
plotSpread(Behavior.GLME.Reffects.PastAct,'xValues',1.15,'distributionColors',[0.52 0.26 0.08],'spreadWidth',0.25);
errorbar([0.85 1.15],Behavior.GLME.Coefficients.Estimate([3 2]),Behavior.GLME.Coefficients.SE([3 2]),'.k','LineWidth',1.5);
xticks([])
ylabel('ME-GLM effect')
legend({'Effect of past outcome (EL)','Effect of past partner''s action (OL)'})

%GLM OL and EL uncertainty
subplot(1,3,2)
hold on
h(1) = plot(0.85:1.85,Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),'-', 'LineWidth', 1.5);
h(2) = plot(1.15:2.15,Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),'-', 'LineWidth', 1.5);
plotSpread([Behavior.GLME_OLunc.Reffects.PastTok_LowOLU Behavior.GLME_OLunc.Reffects.PastTok_HighOLU],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_OLunc.Reffects.PastAct_LowOLU Behavior.GLME_OLunc.Reffects.PastAct_HighOLU],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),Behavior.GLME_OLunc.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),Behavior.GLME_OLunc.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35])
xticks([1 2])
xticklabels({'Low','High'})
xlabel('OL uncertainty')
ylabel('ME-GLM effect')
subplot(1,3,3)
hold on
h(1) = plot(0.85:1.85,Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5);
h(2) = plot(1.15:2.15,Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_ELunc.Reffects.PastTok_LowELU Behavior.GLME_ELunc.Reffects.PastTok_HighELU],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_ELunc.Reffects.PastAct_LowELU Behavior.GLME_ELunc.Reffects.PastAct_HighELU],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),Behavior.GLME_ELunc.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),Behavior.GLME_ELunc.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35])
xticks([1 2])
xticklabels({'Low','High'})
xlabel('EL uncertainty')
ylabel('ME-GLM effect')
legend(h, {'Past outcome (EL effect)','Past action (OL effect)'});


% correlation with actual data
Behavior_actual = load('Behavioral_blockunc.mat');
% general GLME effect
figure;
subplot(1,4,1)
hold on
scatter(Behavior_actual.Behavior.GLME.Reffects.PastAct,Behavior.GLME.Reffects.PastAct,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('PastAct');
refline([1,0]);
subplot(1,4,2)
hold on
scatter(Behavior_actual.Behavior.GLME.Reffects.PastTok,Behavior.GLME.Reffects.PastTok,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('PastTok');
refline([1,0]);
% diff of PastAct between low and high OL unc
diff_OL_actual = Behavior_actual.Behavior.GLME_OLunc.Reffects.PastAct_LowOLU-Behavior_actual.Behavior.GLME_OLunc.Reffects.PastAct_HighOLU;
diff_OL_sim = Behavior.GLME_OLunc.Reffects.PastAct_LowOLU-Behavior.GLME_OLunc.Reffects.PastAct_HighOLU;
subplot(1,4,3)
hold on
scatter(diff_OL_actual,diff_OL_sim,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('OL effect on PastAct');
refline([1,0]);
% diff of PastTok between low and high EL unc
diff_EL_actual = Behavior_actual.Behavior.GLME_ELunc.Reffects.PastTok_LowELU-Behavior_actual.Behavior.GLME_ELunc.Reffects.PastTok_HighELU;
diff_EL_sim = Behavior.GLME_ELunc.Reffects.PastTok_LowELU-Behavior.GLME_ELunc.Reffects.PastTok_HighELU;
subplot(1,4,4)
hold on
scatter(diff_EL_actual,diff_EL_sim,'filled','MarkerFaceAlpha',0.5);
xlabel('actual');
ylabel('simulated');
title('EL effect on PastTok');
refline([1,0]);
