%This script performs all the posterior predictive check analyses for Study 1,
%and generates the plots presented in: Figures 4A, 4C, 4E, 4F, 5C, 6A, 6B,
%S4A, S4B, S4E, S4F, S5A, S5C, S5E, S5G, S8A, S8B

clear all
close all
fs = filesep;

data_dir = ['..' fs '..' fs 'data'];
mod_dir = ['..' fs 'dependencies' fs 'model_functions']; %directory where modelling functions are saved (common to both studies)

addpath(['..' fs 'dependencies']);
addpath(['..' fs 'dependencies' fs 'plotSpread']);
addpath(mod_dir);

% init the randomization screen
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

%load data and format data as needed (in cell array)
data = readtable([data_dir fs 'data_study1.csv']);
subID_list = unique(data.subNb);
n_all = length(subID_list);

%load model fits
load(['model_fitting_outputs' fs 'Recap_model_fitting.mat'])

%load Behavioral variables to compare with model predictions
load('Behavioral_variables.mat')

clr = [248/255 125/255 115/255; 184/255 186/255 65/255; ...
    51/255 198/255 142/255; 34/255 181/255 246/255; ...
    239/255 110/255 253/255];

gsize = [sum(group==1) sum(group==2) sum(group==3) sum(group==4) sum(group==5)];

PPC = struct();

%% Learning curves from model-generated data

%using subjects' best-fitting parameters
learn_sub_Baseline = nan(n_all,8);
learn_sub_EL = nan(n_all,8);
learn_sub_OL = nan(n_all,8);
learn_sub_FixArb = nan(n_all,8);
learn_sub_DynArb = nan(n_all,8);
for s=1:n_all
    params_Baseline = fitRecap.paramRaw.Baseline(s,:);
    params_EL = fitRecap.paramRaw.ExpLearn(s,:);
    params_OL = fitRecap.paramRaw.ObsLearn(s,:);
    params_FixArb = fitRecap.paramRaw.FixArb(s,:);
    params_DynArb = fitRecap.paramRaw.DynArb(s,:);

    subNb = subID_list(s);
    P = table2array(data(data.subNb==subNb,2:end));
    
    %run 1000 iterations to account for stochasticity in choice-generation process
    learn_sim_Baseline = nan(1000,8);
    learn_sim_EL = nan(1000,8);
    learn_sim_OL = nan(1000,8);
    learn_sim_FixArb = nan(1000,8);
    learn_sim_DynArb = nan(1000,8);
    
    for i=1:1000
    
        %predict choice from Baseline model
        P_pred_Baseline = generate_choice_Baseline(params_Baseline, P);
        pred_corr_Baseline = P_pred_Baseline(:,6); % whether choice is correct
        
        %predict choice from EL model
        P_pred_EL = generate_choice_ExpLearn(params_EL, P);
        pred_corr_EL = P_pred_EL(:,6); % whether choice is correct
        
        %predict choice from OL model
        P_pred_OL = generate_choice_ObsLearn(params_OL, P);
        pred_corr_OL = P_pred_OL(:,6); % whether choice is correct
        
        % predict choice from FixArb model
        P_pred_FixArb = generate_choice_FixArb(params_FixArb, P);
        pred_corr_FixArb = P_pred_FixArb(:,6); % whether choice is correct
        
        % predict choice from DynArb model
        P_pred_DynArb = generate_choice_DynArb(params_DynArb, P);
        pred_corr_DynArb = P_pred_DynArb(:,6); % whether choice is correct
       
        % Learning curve
        % identify goal token switch point - count the block break as a break
        switch_idx = [find(diff([0; (P(:,1)*2-2)+P(:,8)]));size(P,1)+1];
        
        learn_mat_Baseline = NaN(length(switch_idx)-1,8);
        learn_mat_EL = NaN(length(switch_idx)-1,8);
        learn_mat_OL = NaN(length(switch_idx)-1,8);
        learn_mat_FixArb = NaN(length(switch_idx)-1,8);
        learn_mat_DynArb = NaN(length(switch_idx)-1,8);

        % only calculate first 8 trials
        for j=1:(length(switch_idx)-1)
            goal_window = [switch_idx(j):switch_idx(j)+7];
            learn_mat_Baseline(j,:) = pred_corr_Baseline(goal_window);
            learn_mat_EL(j,:) = pred_corr_EL(goal_window);
            learn_mat_OL(j,:) = pred_corr_OL(goal_window);
            learn_mat_FixArb(j,:) = pred_corr_FixArb(goal_window);
            learn_mat_DynArb(j,:) = pred_corr_DynArb(goal_window);
        end
        % accuracy for each simulation
        learn_sim_Baseline(i,:) = mean(learn_mat_Baseline,1);
        learn_sim_EL(i,:) = mean(learn_mat_EL,1);
        learn_sim_OL(i,:) = mean(learn_mat_OL,1);
        learn_sim_FixArb(i,:) = mean(learn_mat_FixArb,1);
        learn_sim_DynArb(i,:) = mean(learn_mat_DynArb,1);

    end
    % average for each subject
    learn_sub_Baseline(s,:) = mean(learn_sim_Baseline,1);
    learn_sub_EL(s,:) = mean(learn_sim_EL,1);
    learn_sub_OL(s,:) = mean(learn_sim_OL,1);
    learn_sub_FixArb(s,:) = mean(learn_sim_FixArb,1);
    learn_sub_DynArb(s,:) = mean(learn_sim_DynArb,1);
    
end
PPC.LearningCurves.Baseline = learn_sub_Baseline;
PPC.LearningCurves.ExpLearn = learn_sub_EL;
PPC.LearningCurves.ObsLearn = learn_sub_OL;
PPC.LearningCurves.FixArb = learn_sub_FixArb;
PPC.LearningCurves.DynArb = learn_sub_DynArb;
save('Posterior_predictive_checks.mat','PPC')

%% Plots
%learning curves by groups + model predictions (Figure 5C)
%data
for g=1:5
    mean_learning(:,g) = nanmean(Behavior.Learning(group==g,:))';
    sem_learning(:,g) = nanstd(Behavior.Learning(group==g,:))'/sqrt(gsize(g));
end
%model predictions
mean_learning_pred = [mean(PPC.LearningCurves.Baseline(group==1,:))' mean(PPC.LearningCurves.ExpLearn(group==2,:))' mean(PPC.LearningCurves.ObsLearn(group==3,:))' ...
    mean(PPC.LearningCurves.FixArb(group==4,:))' mean(PPC.LearningCurves.DynArb(group==5,:))'];
%plot
figure; hold
b = plot(mean_learning,'LineWidth',2);
b(1).Color = clr(1,:);
b(2).Color = clr(2,:);
b(3).Color = clr(3,:);
b(4).Color = clr(4,:);
b(5).Color = clr(5,:);
for i=1:5
    shadedErrorBar((1:8), Behavior.Learning(group==i,:), {@mean, @(x) 1*std(x)/sqrt(gsize(i))},'lineProps',{'.','color',clr(i,:)})
end
b = plot(mean_learning_pred,'--','LineWidth',2);
b(1).Color = clr(1,:)-0.1;
b(2).Color = clr(2,:)-0.1;
b(3).Color = clr(3,:)-0.1;
b(4).Color = clr(4,:)-0.1;
b(5).Color = clr(5,:)-0.1;
xlabel('Trial since last reversal')
ylabel('Accuracy')
ylim([0.4 0.85])
legend({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})

%learning curves predicted by each model for shuffled group membership (Figure S8B)
shuffled_groups = sortrows([randperm(n_all)' [repmat(1:5,1,25) randi(5)]'],1);
sgrp = shuffled_groups(:,2);
mean_learning_pred_s = [mean(PPC.LearningCurves.Baseline(sgrp==1,:))' mean(PPC.LearningCurves.ExpLearn(sgrp==2,:))' mean(PPC.LearningCurves.ObsLearn(sgrp==3,:))' ...
    mean(PPC.LearningCurves.FixArb(sgrp==4,:))' mean(PPC.LearningCurves.DynArb(sgrp==5,:))'];
figure; hold
b = plot(mean_learning_pred_s,'--','LineWidth',2);
b(1).Color = clr(1,:);
b(2).Color = clr(2,:);
b(3).Color = clr(3,:);
b(4).Color = clr(4,:);
b(5).Color = clr(5,:);
shadedErrorBar(1:8, PPC.LearningCurves.Baseline(sgrp==1,:), {@mean, @(x) 1*std(x)/sqrt(sum(sgrp==1))},'lineProps',{'.','color',clr(1,:)})
shadedErrorBar(1:8, PPC.LearningCurves.ExpLearn(sgrp==2,:), {@mean, @(x) 1*std(x)/sqrt(sum(sgrp==2))},'lineProps',{'.','color',clr(2,:)})
shadedErrorBar(1:8, PPC.LearningCurves.ObsLearn(sgrp==3,:), {@mean, @(x) 1*std(x)/sqrt(sum(sgrp==3))},'lineProps',{'.','color',clr(3,:)})
shadedErrorBar(1:8, PPC.LearningCurves.FixArb(sgrp==4,:), {@mean, @(x) 1*std(x)/sqrt(sum(sgrp==4))},'lineProps',{'.','color',clr(4,:)})
shadedErrorBar(1:8, PPC.LearningCurves.DynArb(sgrp==5,:), {@mean, @(x) 1*std(x)/sqrt(sum(sgrp==5))},'lineProps',{'.','color',clr(5,:)})
xlabel('Trial since last reversal')
ylabel('Accuracy')
ylim([0.4 0.85])
legend({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
title('Model prediction by randomly shuffled group')


%% Learning curve using simulated parameters similar across models (Figure S8A)
lr_arb = 0.05:0.1:0.95;
lr_arb_r = -log((1-lr_arb)./lr_arb);
i=1;
for a=1:10
    for b=1:10
        lr_arb_all(i,:) = [lr_arb_r(a) lr_arb_r(b)];
        i = i+1;
    end
end
w = sortrows([randperm(100)' linspace(-1,1)']);
%mag_boost = sortrows([randperm(100)' linspace(-1,5)']);
%keep beta and mag_boost constant

%select 1 trial list of each type from P matrix (from study 1 data summary, 
%16 different trial lists in total)
ind_trial_list = [10;27;8;6;18;16;5;2;12;4;25;29;19;9;20;21];

learn_sub_EL = nan(100,8);
learn_sub_OL = nan(100,8);
learn_sub_FixArb = nan(100,8);
learn_sub_DynArb = nan(100,8);
recap_perf = [1./(1+exp(-lr_arb_all)) 1./(1+exp(-w(:,2))) nan(100,4)];
for s=1:100
    s
    params_EL = [2 lr_arb_all(s,1) 1];
    params_OL = [2 lr_arb_all(s,2)];
    params_FixArb = [2 2 lr_arb_all(s,1) lr_arb_all(s,2) w(s,2) 1];
    params_DynArb = [2 2 lr_arb_all(s,1) lr_arb_all(s,2) 1 w(s,2)];
    
    %select one of the 16 trial lists at random
    tr_list = ind_trial_list(randi(16));
    P = table2array(data(data.subNb == tr_list,2:end));
    
    %run 1000 iterations to account for stochasticity in choice-generation process
    learn_sim_Baseline = nan(1000,8);
    learn_sim_EL = nan(1000,8);
    learn_sim_OL = nan(1000,8);
    learn_sim_FixArb = nan(1000,8);
    learn_sim_DynArb = nan(1000,8);
    model_perf = nan(1000,4);
    for i=1:1000
  
        %predict choice from EL model
        P_pred_EL = generate_choice_ExpLearn(params_EL, P);
        pred_corr_EL = P_pred_EL(:,6); % whether choice is correct
        
        %predict choice from OL model
        P_pred_OL = generate_choice_ObsLearn(params_OL, P);
        pred_corr_OL = P_pred_OL(:,6); % whether choice is correct
        
        % predict choice from FixArb model
        P_pred_FixArb = generate_choice_FixArb(params_FixArb, P);
        pred_corr_FixArb = P_pred_FixArb(:,6); % whether choice is correct
        
        % predict choice from DynArb model
        P_pred_DynArb = generate_choice_DynArb(params_DynArb, P);
        pred_corr_DynArb = P_pred_DynArb(:,6); % whether choice is correct
        
        % accuracy recap
        model_perf(i,:) = [nanmean(P_pred_EL(:,6)) nanmean(P_pred_OL(:,6)) ...
            nanmean(P_pred_FixArb(:,6)) nanmean(P_pred_DynArb(:,6))];
        
        % Learning curve
        % identify goal token switch point - count the block break as a break
        switch_idx = [find(diff([0; (P(:,1)*2-2)+P(:,8)]));size(P,1)+1];
        
        learn_mat_EL = NaN(length(switch_idx)-1,8);
        learn_mat_OL = NaN(length(switch_idx)-1,8);
        learn_mat_FixArb = NaN(length(switch_idx)-1,8);
        learn_mat_DynArb = NaN(length(switch_idx)-1,8);

        % only calculate first 8 trials
        for j=1:(length(switch_idx)-1)
            goal_window = [switch_idx(j):switch_idx(j)+7];
            learn_mat_EL(j,:) = pred_corr_EL(goal_window);
            learn_mat_OL(j,:) = pred_corr_OL(goal_window);
            learn_mat_FixArb(j,:) = pred_corr_FixArb(goal_window);
            learn_mat_DynArb(j,:) = pred_corr_DynArb(goal_window);
        end
        % accuracy for each simulation
        learn_sim_EL(i,:) = mean(learn_mat_EL,1);
        learn_sim_OL(i,:) = mean(learn_mat_OL,1);
        learn_sim_FixArb(i,:) = mean(learn_mat_FixArb,1);
        learn_sim_DynArb(i,:) = mean(learn_mat_DynArb,1);

    end
    % average for each subject
    learn_sub_EL(s,:) = mean(learn_sim_EL,1);
    learn_sub_OL(s,:) = mean(learn_sim_OL,1);
    learn_sub_FixArb(s,:) = mean(learn_sim_FixArb,1);
    learn_sub_DynArb(s,:) = mean(learn_sim_DynArb,1);
    recap_perf(s,4:7) = mean(model_perf,1);
    
end
PPC.LearningCurvesSimParams.ExpLearn = learn_sub_EL;
PPC.LearningCurvesSimParams.ObsLearn = learn_sub_OL;
PPC.LearningCurvesSimParams.FixArb = learn_sub_FixArb;
PPC.LearningCurvesSimParams.DynArb = learn_sub_DynArb;
PPC.LearningCurvesSimParams.Recap = recap_perf;
save('Posterior_predictive_checks.mat','PPC')

%plot learning curve, simulated parameters (100 sets) (Figure S8A)
mean_learning_pred = [mean(learn_sub_EL)' mean(learn_sub_OL)' mean(learn_sub_FixArb)' mean(learn_sub_DynArb)'];
figure; hold
b = plot(mean_learning_pred,'--','LineWidth',2);
b(1).Color = clr(2,:);
b(2).Color = clr(3,:);
b(3).Color = clr(4,:);
b(4).Color = clr(5,:);
shadedErrorBar(1:8, learn_sub_EL, {@mean, @(x) 1*std(x)/sqrt(100)},'lineProps',{'.','color',clr(2,:)})
shadedErrorBar(1:8, learn_sub_OL, {@mean, @(x) 1*std(x)/sqrt(100)},'lineProps',{'.','color',clr(3,:)})
shadedErrorBar(1:8, learn_sub_FixArb, {@mean, @(x) 1*std(x)/sqrt(100)},'lineProps',{'.','color',clr(4,:)})
shadedErrorBar(1:8, learn_sub_DynArb, {@mean, @(x) 1*std(x)/sqrt(100)},'lineProps',{'.','color',clr(5,:)})
xlabel('Trial since last reversal')
ylabel('Accuracy')
ylim([0.4 0.85])
legend({'ExpLearn','ObsLearn','FixArb','DynArb'})
title('Model prediction (simulated parameters)')

%fix vs dyn arb model performance 
recap_perf_w = sort(recap_perf(:,[3 4:7]),1);
figure; hold
plot(recap_perf_w(:,1),recap_perf_w(:,4),'--','Color',clr(4,:),'LineWidth',2)
plot(recap_perf_w(:,1),recap_perf_w(:,5),'--','Color',clr(5,:),'LineWidth',2)

perf_diff = recap_perf(:,7) - recap_perf(:,6);
recap_mat = nan(10,10);
i=1;
for r = 1:10
    for c = 1:10 
        recap_mat(r,c) = perf_diff(i);
        i=i+1;
    end
end
figure;
h = heatmap(unique(recap_perf(:,1)), unique(recap_perf(:,2)), recap_mat);
h.Title = 'Performance difference (dynamic vs fixed arbitration model)';
h.XLabel = 'EL learning rate';
h.YLabel = 'OL learning rate';
h.Colormap = redwhiteblue(0,1);
h.ColorLimits = [0 0.12];

perf_diff = recap_perf(:,5) - recap_perf(:,4);
recap_mat = nan(10,10);
i=1;
for r = 1:10
    for c = 1:10 
        recap_mat(r,c) = perf_diff(i);
        i=i+1;
    end
end
figure;
h = heatmap(unique(recap_perf(:,1)), unique(recap_perf(:,2)), recap_mat);
h.Title = 'Performance difference (OL vs EL model)';
h.XLabel = 'EL learning rate';
h.YLabel = 'OL learning rate';
h.Colormap = redwhiteblue(-1,1);
h.ColorLimits = [-0.2 0.2];

%% Check that OL vs EL models predict corresponding behavioral signature (Figure S4)
ntr = 160;
nt_diff_from_OL_sub = nan(n_all,1);
prop_OL_ch_from_OL_sub = nan(n_all,1);
nt_diff_from_EL_sub = nan(n_all,1);
prop_OL_ch_from_EL_sub = nan(n_all,1);
for s=1:n_all
    s
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
PPC.OLch.ExpLearnPredCh = prop_OL_ch_from_EL_sub;
PPC.OLch.ObsLearnPredCh = prop_OL_ch_from_OL_sub;
PPC.OLch.ExpLearnNTdiff = nt_diff_from_EL_sub;
PPC.OLch.ObsLearnNTdiff = nt_diff_from_OL_sub;
save('Posterior_predictive_checks.mat','PPC')

%plot histograms of model predictions and data (Figure S4A-B)
figure;
subplot(2,1,1); hold on
histogram(PPC.OLch.ExpLearnPredCh,20,'FaceAlpha',0.6);
histogram(PPC.OLch.ObsLearnPredCh,15,'FaceAlpha',0.6);
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

%plot correlations between data and model predictions (Figure S4E-F)
iEL = Behavior.Prop_OL_ch(:,1)<0.5;
iOL = Behavior.Prop_OL_ch(:,1)>0.5;
figure;
subplot(1,2,1); hold on
plot(Behavior.Prop_OL_ch(iEL,1),PPC.OLch.ExpLearnPredCh(iEL),'.','Color','#0072BD'); lsline()
plot(Behavior.Prop_OL_ch(iEL,1),PPC.OLch.ObsLearnPredCh(iEL),'.','Color','#D95319'); lsline()
xlim([0.2 0.5]); ylim([0.2 0.8])
title({'Experiential learners'; '(OL choice prop. < 0.5)'})
xlabel('Data')
ylabel('Model predictions')
legend({'EL model','','OL model',''})
subplot(1,2,2); hold on
plot(Behavior.Prop_OL_ch(iOL,1),PPC.OLch.ExpLearnPredCh(iOL),'.','Color','#0072BD'); lsline()
plot(Behavior.Prop_OL_ch(iOL,1),PPC.OLch.ObsLearnPredCh(iOL),'.','Color','#D95319'); lsline()
xlim([0.5 0.8]); ylim([0.2 0.8])
title({'Observational learners'; '(OL choice prop. > 0.5)'})
xlabel('Data')

%print correlations
[r,p] = corr(Behavior.Prop_OL_ch(iEL,1),PPC.OLch.ExpLearnPredCh(iEL))
[r,p] = corr(Behavior.Prop_OL_ch(iEL,1),PPC.OLch.ObsLearnPredCh(iEL))
[r,p] = corr(Behavior.Prop_OL_ch(iOL,1),PPC.OLch.ExpLearnPredCh(iOL))
[r,p] = corr(Behavior.Prop_OL_ch(iOL,1),PPC.OLch.ObsLearnPredCh(iOL))

%% Run "main effect" glm on data generated by OL, EL, and DynArb
sim_mod = {'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'};
gc_mod_list = {@generate_choice_Baseline; @generate_choice_ExpLearn; @generate_choice_ObsLearn; ...
    @generate_choice_FixArb; @generate_choice_DynArb};

pred_EL_effect_subs = nan(n_all,5);
pred_OL_effect_subs = nan(n_all,5);
pred_EL_effect_fe = nan(2,5);
pred_OL_effect_fe = nan(2,5);

for m=1:5
    params_all = fitRecap.paramRaw.(sim_mod{m}); 
    glm_mat_allsubs = [];
    for sub=1:n_all

        %load subject data
        subNb = Behavior.subID_list(sub);
        P = table2array(data(data.subNb==subNb,2:end));

        %generate data
        params_sub = params_all(sub,:);
        P_pred = gc_mod_list{m}(params_sub, P);

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

    end

    %run glm
    glm_tbl = array2table(glm_mat_allsubs, 'VariableNames',{'subNb','PastAct','PastTok','choice'});
    glme = fitglme(glm_tbl,'choice ~ 1 + PastAct + PastTok + (1 + PastAct + PastTok|subNb)', ...
        'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    [~,~,stats] = randomEffects(glme);
    Meffects = dataset2table(glme.Coefficients);
    Reffects = dataset2table(stats);
    U = unstack(Reffects(:,2:4),'Estimate','Name');
    U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
    pred_OL_effect_subs(:,m) = U.PastAct + Meffects.Estimate(2);
    pred_EL_effect_subs(:,m) = U.PastTok + Meffects.Estimate(3);
    pred_OL_effect_fe(:,m) = [Meffects.Estimate(2); Meffects.SE(2)];
    pred_EL_effect_fe(:,m) = [Meffects.Estimate(3); Meffects.SE(3)];
        
end
PPC.GLME_main.PastTok = pred_EL_effect_subs;
PPC.GLME_main.PastAct = pred_OL_effect_subs;
PPC.GLME_main.PastTok_FE = pred_EL_effect_fe;
PPC.GLME_main.PastAct_FE = pred_OL_effect_fe;
save('Posterior_predictive_checks.mat','PPC')

%% Plots for Figure 4
%bar plot of mean actual and predicted GLM effects (Figure 4A)
figure;
subplot(1,2,1); hold on
boxchart([Behavior.GLME.Reffects.PastTok PPC.GLME_main.PastTok],'MarkerStyle','.')
plot([0 0 0 0 0 0 0], 'k--')
xticklabels({'Data','Baseline model','EL model','OL model','FixArb model','DynArb model'})
xtickangle(30)
title('Effect of past outcome (EL)')
subplot(1,2,2); hold on
boxchart([Behavior.GLME.Reffects.PastAct PPC.GLME_main.PastAct],'MarkerStyle','.',...
    'MarkerColor',"#D95319",'BoxFaceColor',"#D95319")
plot([0 0 0 0 0 0 0], 'k--')
xticklabels({'Data','Baseline model','EL model','OL model','FixArb model','DynArb model'})
xtickangle(30)
title('Effect of past partner''s action (OL)')

[~,p,~,s] = ttest(PPC.GLME_main.PastTok(:,3),Behavior.GLME.Reffects.PastTok) %EL effect for OL model vs behavior
[~,p,~,s] = ttest(PPC.GLME_main.PastTok(:,3),PPC.GLME_main.PastTok(:,2)) %EL effect for OL model vs EL model
[~,p,~,s] = ttest(PPC.GLME_main.PastTok(:,3),PPC.GLME_main.PastTok(:,5)) %EL effect for OL model vs DynArb model
[~,p,~,s] = ttest(PPC.GLME_main.PastAct(:,2),Behavior.GLME.Reffects.PastAct) %OL effect for EL model vs behavior
[~,p,~,s] = ttest(PPC.GLME_main.PastAct(:,2),PPC.GLME_main.PastAct(:,3)) %OL effect for EL model vs OL model
[~,p,~,s] = ttest(PPC.GLME_main.PastAct(:,2),PPC.GLME_main.PastAct(:,5)) %OL effect for EL model vs DynArb model

%scatter plot to show PPC across participants (Figure 4C)
EL_col = [0 0.4470 0.7410];
OL_col = [0.8500 0.3250 0.0980];
figure;
subplot(2,3,1); hold on
plot(Behavior.GLME.Reffects.PastTok,PPC.GLME_main.PastTok(:,2),'.','MarkerEdgeColor',EL_col)
lsline()
xlabel('EL effect - data'); ylabel('EL model')
ylim([-0.5 2])
title(['R=' num2str(corr(Behavior.GLME.Reffects.PastTok,PPC.GLME_main.PastTok(:,2)),3)])
subplot(2,3,2); hold on
plot(Behavior.GLME.Reffects.PastTok,PPC.GLME_main.PastTok(:,3),'.','MarkerEdgeColor',EL_col)
lsline()
xlabel('EL effect - data'); ylabel('OL model')
ylim([-0.5 2])
title(['R=' num2str(corr(Behavior.GLME.Reffects.PastTok,PPC.GLME_main.PastTok(:,3)),3)])
subplot(2,3,3); hold on
plot(Behavior.GLME.Reffects.PastTok,PPC.GLME_main.PastTok(:,5),'.','MarkerEdgeColor',EL_col)
lsline()
xlabel('EL effect - data'); ylabel('DynArb model')
ylim([-0.5 2])
title(['R=' num2str(corr(Behavior.GLME.Reffects.PastTok,PPC.GLME_main.PastTok(:,5)),3)])
subplot(2,3,4); hold on
plot(Behavior.GLME.Reffects.PastAct,PPC.GLME_main.PastAct(:,2),'.','MarkerEdgeColor',OL_col)
h = lsline(); h.Color = OL_col;
xlabel('OL effect - data'); ylabel('EL model')
ylim([-0.1 0.7])
title(['R=' num2str(corr(Behavior.GLME.Reffects.PastAct,PPC.GLME_main.PastAct(:,2)),3)])
subplot(2,3,5); hold on
plot(Behavior.GLME.Reffects.PastAct,PPC.GLME_main.PastAct(:,3),'.','MarkerEdgeColor',OL_col)
h = lsline(); h.Color = OL_col;
xlabel('OL effect - data'); ylabel('OL model')
ylim([-0.1 0.7])
title(['R=' num2str(corr(Behavior.GLME.Reffects.PastAct,PPC.GLME_main.PastAct(:,3)),3)])
subplot(2,3,6); hold on
plot(Behavior.GLME.Reffects.PastAct,PPC.GLME_main.PastAct(:,5),'.','MarkerEdgeColor',OL_col)
h = lsline(); h.Color = OL_col;
xlabel('OL effect - data'); ylabel('DynArb model')
ylim([-0.1 0.7])
title(['R=' num2str(corr(Behavior.GLME.Reffects.PastAct,PPC.GLME_main.PastAct(:,5)),3)])


%% Breakdown main glm effects by group
%data
mean_glm_data = nan(2,5); sem_glm_data = nan(2,5);
for g=1:5    
    mean_glm_data(:,g) = [nanmean(Behavior.GLME.Reffects.PastTok(group==g)); ...
        nanmean(Behavior.GLME.Reffects.PastAct(group==g))];
    sem_glm_data(:,g) = [nanstd(Behavior.GLME.Reffects.PastTok(group==g))/sqrt(gsize(g)); ...
        nanstd(Behavior.GLME.Reffects.PastAct(group==g))/sqrt(gsize(g))];
end

%using original group definition and best-fitting model for each group
pred_effects = nan(n_all,2);
for s=1:n_all
    pred_effects(s,:)=[PPC.GLME_main.PastTok(s,group(s)) PPC.GLME_main.PastAct(s,group(s))];
end

mean_glm = [mean(pred_effects(group==1,:));mean(pred_effects(group==2,:));mean(pred_effects(group==3,:));...
    mean(pred_effects(group==4,:));mean(pred_effects(group==5,:))];
sem_glm = [std(pred_effects(group==1,:))/sqrt(sum(group==1));std(pred_effects(group==2,:))/sqrt(sum(group==2));...
    std(pred_effects(group==3,:))/sqrt(sum(group==3));...
    std(pred_effects(group==4,:))/sqrt(sum(group==4));std(pred_effects(group==5,:))/sqrt(sum(group==5))];

figure;
subplot(1,2,1); hold on
bar(1:5,mean_glm(:,1),'EdgeColor','k','LineWidth',1);
plotSpread(pred_effects(:,1),'distributionIdx',group,'xValues',(1:5),'distributionColors',[0 0.32 0.47]);
errorbar((1:5),mean_glm(:,1),sem_glm(:,1),'.k','LineWidth',1.5);
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
ylim([-0.5 2])
xlabel('Group')
ylabel('ME-GLM effect')
title(["Model prediction by group","EL effect"])
subplot(1,2,2); hold on
bar(1:5,mean_glm(:,2),'FaceColor','#D95319','EdgeColor','k','LineWidth',1);
plotSpread(pred_effects(:,2),'distributionIdx',group,'xValues',(1:5),'distributionColors',[0.52 0.26 0.08]);
errorbar((1:5),mean_glm(:,2),sem_glm(:,2),'.k','LineWidth',1.5);
ylim([-0.3 1])
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('ME-GLM effect')
title(["Model prediction by group","OL effect"])

%randomly shuffling group membership
shuffled_groups = sortrows([randperm(n_all)' [repmat(1:5,1,25) randi(5)]'],1);
sgrp = shuffled_groups(:,2);

pred_effects_sh = nan(n_all,2);
for s=1:n_all
    pred_effects_sh(s,:)=[PPC.GLME_main.PastTok(s,sgrp(s)) PPC.GLME_main.PastAct(s,sgrp(s))];
end

mean_glm_sh = [mean(pred_effects_sh(sgrp==1,:));mean(pred_effects_sh(sgrp==2,:));mean(pred_effects_sh(sgrp==3,:));...
    mean(pred_effects_sh(sgrp==4,:));mean(pred_effects_sh(sgrp==5,:))];
sem_glm_sh = [std(pred_effects_sh(sgrp==1,:))/sqrt(sum(sgrp==1));std(pred_effects_sh(sgrp==2,:))/sqrt(sum(sgrp==2));...
    std(pred_effects_sh(sgrp==3,:))/sqrt(sum(sgrp==3));...
    std(pred_effects_sh(sgrp==4,:))/sqrt(sum(sgrp==4));std(pred_effects_sh(sgrp==5,:))/sqrt(sum(sgrp==5))];

figure;
subplot(1,2,1); hold on
bar(1:5,mean_glm_sh(:,1),'EdgeColor','k','LineWidth',1);
plotSpread(pred_effects_sh(:,1),'distributionIdx',sgrp,'xValues',(1:5),'distributionColors',[0 0.32 0.47]);
errorbar((1:5),mean_glm_sh(:,1),sem_glm_sh(:,1),'.k','LineWidth',1.5);
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
ylim([-0.5 2])
xlabel('Group')
ylabel('ME-GLM effect')
title(["Model prediction by randomly shuffled group","EL effect"])
subplot(1,2,2); hold on
bar(1:5,mean_glm_sh(:,2),'FaceColor','#D95319','EdgeColor','k','LineWidth',1);
plotSpread(pred_effects_sh(:,2),'distributionIdx',sgrp,'xValues',(1:5),'distributionColors',[0.52 0.26 0.08]);
errorbar((1:5),mean_glm_sh(:,2),sem_glm_sh(:,2),'.k','LineWidth',1.5);
ylim([-0.3 1])
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('ME-GLM effect')
title(["Model prediction by randomly shuffled group","OL effect"])

%% Plot Figure 6A-B
figure;
subplot(1,2,1); hold on
bar(1:3:13,mean_glm_data(1,:),0.25,'EdgeColor','k','FaceColor',[79 184 255]/255,'LineWidth',1);
bar((1:3:13)+0.8, mean_glm(:,1),0.25,'EdgeColor','k','FaceColor',[0 114 189]/255,'LineWidth',1);
bar((1:3:13)+1.6,mean_glm_sh(:,1),0.25,'EdgeColor','k','FaceColor',[0.6 0.6 0.6],'LineWidth',1);
plotSpread(Behavior.GLME.Reffects.PastTok,'distributionIdx',group,'xValues',(1:3:13),'distributionColors',[21 160 255]/255);
errorbar((1:3:13),mean_glm_data(1,:),sem_glm_data(1,:),'.k','LineWidth',1,'CapSize',3);
errorbar((1:3:13)+0.8,mean_glm(:,1),sem_glm(:,1),'.k','LineWidth',1,'CapSize',3);
errorbar((1:3:13)+1.6,mean_glm_sh(:,1),sem_glm_sh(:,1),'.k','LineWidth',1,'CapSize',3);
xticks(2:3:14); xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'}); xtickangle(30); xlabel('Group')
ylim([-0.5 2]); ylabel('ME-GLM effect')
legend({'data','model by group','model by shuffled group'})
title('Effect of past outcome (EL)')
subplot(1,2,2); hold on
bar(1:3:13,mean_glm_data(2,:),0.25,'EdgeColor','k','FaceColor',[246 146 100]/255,'LineWidth',1);
bar((1:3:13)+0.8, mean_glm(:,2),0.25,'EdgeColor','k','FaceColor',[176 60 8]/255,'LineWidth',1);
bar((1:3:13)+1.6, mean_glm_sh(:,2),0.25,'EdgeColor','k','FaceColor',[0.6 0.6 0.6],'LineWidth',1);
plotSpread(Behavior.GLME.Reffects.PastAct,'distributionIdx',group,'xValues',(1:3:13),'distributionColors',[246 115 56]/255);
errorbar((1:3:13),mean_glm_data(2,:),sem_glm_data(2,:),'.k','LineWidth',1,'CapSize',3);
errorbar((1:3:13)+0.8,mean_glm(:,2),sem_glm(:,2),'.k','LineWidth',1,'CapSize',3);
errorbar((1:3:13)+1.6,mean_glm_sh(:,2),sem_glm_sh(:,2),'.k','LineWidth',1,'CapSize',3);
ylim([-0.3 1])
xticks(2:3:14); xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'}); xtickangle(30); xlabel('Group')
legend({'data','model by group','model by shuffled group'})
title('Effect of past partner''s action (OL)')




%% Run mixed effect glms with uncertainty effects on data generated by OL, EL, and DynArb
%use both definition of uncertainty: uncertainty trials and uncertainty condition (design blocks)
sim_mod = {'ExpLearn','ObsLearn','DynArb'};
gc_mod_list = {@generate_choice_ExpLearn; @generate_choice_ObsLearn; @generate_choice_DynArb};

ntr = 160;
for m=1:3
    params_all = fitRecap.paramRaw.(sim_mod{m}); 
    glm_mat_OLunc = [];
    glm_mat_ELunc = [];
    glm_mat_Mag = [];
    glm_mat_OLunc_des = [];
    glm_mat_ELunc_des = [];
    
    for sub=1:n_all

        %load subject data
        subNb = Behavior.subID_list(sub);
        P = table2array(data(data.subNb==subNb,2:end));
        
        %generate data
        params_sub = params_all(sub,:);
        P_pred = gc_mod_list{m}(params_sub, P);

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
        
        %define OL and EL uncertainty conditions from design for comparison
        OL_lu_design = P(:,5) > 0.7;
        EL_lu_design = P(:,7) > 0.7 | P(:,7) < 0.3;

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

        %build glm matrix to look at interactions with uncertainty in mixed-effect GLM
        glm_mat_OLunc = [glm_mat_OLunc; ones(sum(good),1)*subNb ...
            [zscore(past_act(good & OL_lu)); zeros(sum(good & ~OL_lu),1)] ...
            [zeros(sum(good & OL_lu),1); zscore(past_act(good & ~OL_lu))] ...
            [zscore(past_tok_comb(good & OL_lu)); zeros(sum(good & ~OL_lu),1)] ...
            [zeros(sum(good & OL_lu),1); zscore(past_tok_comb(good & ~OL_lu))] ...
            [P_pred(good & OL_lu,4);P_pred(good & ~OL_lu,4)]];
        
        glm_mat_OLunc_des = [glm_mat_OLunc_des; ones(sum(good),1)*subNb ...
            [zscore(past_act(good & OL_lu_design)); zeros(sum(good & ~OL_lu_design),1)] ...
            [zeros(sum(good & OL_lu_design),1); zscore(past_act(good & ~OL_lu_design))] ...
            [zscore(past_tok_comb(good & OL_lu_design)); zeros(sum(good & ~OL_lu_design),1)] ...
            [zeros(sum(good & OL_lu_design),1); zscore(past_tok_comb(good & ~OL_lu_design))] ...
            [P_pred(good & OL_lu_design,4);P_pred(good & ~OL_lu_design,4)]];

        glm_mat_ELunc = [glm_mat_ELunc; ones(sum(good),1)*subNb ...
            [zscore(past_act(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
            [zeros(sum(good & EL_lu),1); zscore(past_act(good & ~EL_lu))] ...
            [zscore(past_tok_comb(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
            [zeros(sum(good & EL_lu),1); zscore(past_tok_comb(good & ~EL_lu))] ...
            [P_pred(good & EL_lu,4);P_pred(good & ~EL_lu,4)]];
        
        glm_mat_ELunc_des = [glm_mat_ELunc_des; ones(sum(good),1)*subNb ...
            [zscore(past_act(good & EL_lu_design)); zeros(sum(good & ~EL_lu_design),1)] ...
            [zeros(sum(good & EL_lu_design),1); zscore(past_act(good & ~EL_lu_design))] ...
            [zscore(past_tok_comb(good & EL_lu_design)); zeros(sum(good & ~EL_lu_design),1)] ...
            [zeros(sum(good & EL_lu_design),1); zscore(past_tok_comb(good & ~EL_lu_design))] ...
            [P_pred(good & EL_lu_design,4);P_pred(good & ~EL_lu_design,4)]];

        glm_mat_Mag = [glm_mat_Mag; ones(sum(good),1)*subNb ...
            [zscore(past_act(good & RM_h)); zeros(sum(good & ~RM_h),1)] ...
            [zeros(sum(good & RM_h),1); zscore(past_act(good & ~RM_h))] ...
            [zscore(past_tok_comb(good & RM_h)); zeros(sum(good & ~RM_h),1)] ...
            [zeros(sum(good & RM_h),1); zscore(past_tok_comb(good & ~RM_h))] ...
            [P_pred(good & RM_h,4);P_pred(good & ~RM_h,4)]];

    end

    %separate GLMs for each uncertainty condition
    %OL uncertainty by uncertainty trial type
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
    PPC.GLME_OLunc.PastAct_LowOLU(:,m) = U.PastAct_LowOLU + Meffects.Estimate(2);
    PPC.GLME_OLunc.PastAct_HighOLU(:,m) = U.PastAct_HighOLU + Meffects.Estimate(3);
    PPC.GLME_OLunc.PastTok_LowOLU(:,m) = U.PastTok_LowOLU + Meffects.Estimate(4);
    PPC.GLME_OLunc.PastTok_HighOLU(:,m) = U.PastTok_HighOLU + Meffects.Estimate(5);

    %OL uncertainty by uncertainty condition from the design
    glm_OLunc_des_tbl = array2table(glm_mat_OLunc_des, 'VariableNames',...
        {'subNb','PastAct_LowOLU','PastAct_HighOLU','PastTok_LowOLU','PastTok_HighOLU','choice'});
    glme_OLunc_des = fitglme(glm_OLunc_des_tbl,...
        'choice ~ 1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU + (1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU|subNb)', ...
        'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    %disp(glme_OLunc)
    [~,~,stats] = randomEffects(glme_OLunc_des);
    Meffects = dataset2table(glme_OLunc_des.Coefficients);
    Reffects = dataset2table(stats);
    U = unstack(Reffects(:,2:4),'Estimate','Name');
    PPC.GLME_OLunc_des.PastAct_LowOLU(:,m) = U.PastAct_LowOLU + Meffects.Estimate(2);
    PPC.GLME_OLunc_des.PastAct_HighOLU(:,m) = U.PastAct_HighOLU + Meffects.Estimate(3);
    PPC.GLME_OLunc_des.PastTok_LowOLU(:,m) = U.PastTok_LowOLU + Meffects.Estimate(4);
    PPC.GLME_OLunc_des.PastTok_HighOLU(:,m) = U.PastTok_HighOLU + Meffects.Estimate(5);

    %EL uncertainty by uncertainty trial type
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
    PPC.GLME_ELunc.PastAct_LowELU(:,m) = U.PastAct_LowELU + Meffects.Estimate(2);
    PPC.GLME_ELunc.PastAct_HighELU(:,m) = U.PastAct_HighELU + Meffects.Estimate(3);
    PPC.GLME_ELunc.PastTok_LowELU(:,m) = U.PastTok_LowELU + Meffects.Estimate(4);
    PPC.GLME_ELunc.PastTok_HighELU(:,m) = U.PastTok_HighELU + Meffects.Estimate(5);

    %EL uncertainty by uncertainty condition from the design
    glm_ELunc_des_tbl = array2table(glm_mat_ELunc_des, 'VariableNames',...
        {'subNb','PastAct_LowELU','PastAct_HighELU','PastTok_LowELU','PastTok_HighELU','choice'});
    glme_ELunc_des = fitglme(glm_ELunc_des_tbl,...
        'choice ~ 1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU + (1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU|subNb)', ...
        'Distribution','Binomial','Link','logit','FitMethod','Laplace');
    %disp(glme_ELunc)
    [~,~,stats] = randomEffects(glme_ELunc_des);
    Meffects = dataset2table(glme_ELunc_des.Coefficients);
    Reffects = dataset2table(stats);
    U = unstack(Reffects(:,2:4),'Estimate','Name');
    PPC.GLME_ELunc_des.PastAct_LowELU(:,m) = U.PastAct_LowELU + Meffects.Estimate(2);
    PPC.GLME_ELunc_des.PastAct_HighELU(:,m) = U.PastAct_HighELU + Meffects.Estimate(3);
    PPC.GLME_ELunc_des.PastTok_LowELU(:,m) = U.PastTok_LowELU + Meffects.Estimate(4);
    PPC.GLME_ELunc_des.PastTok_HighELU(:,m) = U.PastTok_HighELU + Meffects.Estimate(5);

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
    PPC.GLME_Mag.PastAct_HighMag(:,m) = U.PastAct_HighMag + Meffects.Estimate(2);
    PPC.GLME_Mag.PastAct_LowMag(:,m) = U.PastAct_LowMag + Meffects.Estimate(3);
    PPC.GLME_Mag.PastTok_HighMag(:,m) = U.PastTok_HighMag + Meffects.Estimate(4);
    PPC.GLME_Mag.PastTok_LowMag(:,m) = U.PastTok_LowMag + Meffects.Estimate(5);

end
save('Posterior_predictive_checks.mat','PPC')


%GLM effect of OL uncertainty (trial definition)
recap_pred_EL = [PPC.GLME_OLunc.PastTok_LowOLU(:,1) PPC.GLME_OLunc.PastTok_HighOLU(:,1) ...
    PPC.GLME_OLunc.PastAct_LowOLU(:,1) PPC.GLME_OLunc.PastAct_HighOLU(:,1)];
recap_pred_OL = [PPC.GLME_OLunc.PastTok_LowOLU(:,2) PPC.GLME_OLunc.PastTok_HighOLU(:,2) ...
    PPC.GLME_OLunc.PastAct_LowOLU(:,2) PPC.GLME_OLunc.PastAct_HighOLU(:,2)];
recap_pred_Arb = [PPC.GLME_OLunc.PastTok_LowOLU(:,3) PPC.GLME_OLunc.PastTok_HighOLU(:,3) ...
    PPC.GLME_OLunc.PastAct_LowOLU(:,3) PPC.GLME_OLunc.PastAct_HighOLU(:,3)];
figure;
subplot(1,4,1); hold on
plot(0.85:1.85,Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5);
plot(1.15:2.15,Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_OLunc.Reffects.PastTok_LowOLU Behavior.GLME_OLunc.Reffects.PastTok_HighOLU],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_OLunc.Reffects.PastAct_LowOLU Behavior.GLME_OLunc.Reffects.PastAct_HighOLU],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),Behavior.GLME_OLunc.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),Behavior.GLME_OLunc.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1.5 3])
xlabel('OL uncertainty')
ylabel('ME-GLM effect')
title('Data')
subplot(1,4,2); hold on
plot(0.85:1.85,mean(recap_pred_EL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_EL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_EL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_EL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_EL(:,1:2)),std(recap_pred_EL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_EL(:,3:4)),std(recap_pred_EL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1.5 3])
xlabel('OL uncertainty')
title('EL model')
subplot(1,4,3); hold on
plot(0.85:1.85,mean(recap_pred_OL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_OL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_OL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_OL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_OL(:,1:2)),std(recap_pred_OL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_OL(:,3:4)),std(recap_pred_OL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1.5 3])
xlabel('OL uncertainty')
title('OL model')
subplot(1,4,4); hold on
plot(0.85:1.85,mean(recap_pred_Arb(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_Arb(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_Arb(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_Arb(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_Arb(:,1:2)),std(recap_pred_Arb(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_Arb(:,3:4)),std(recap_pred_Arb(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1.5 3])
xlabel('OL uncertainty')
title('DynArb model')
legend({'Past outcome (EL effect)','Past action (OL effect)'})

% Simplified plot (Figure 4E)
figure;
plot([0.8 1.2],Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5, 'Color',[21 160 255]/255); hold on
plot([1.8 2.2],Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5, 'Color',[246 115 56]/255);
plot([0.85 1.25],mean(recap_pred_Arb(:,1:2)),'--', 'LineWidth', 1.5, 'Color',[0 114 189]/255);
plot([1.85 2.25],mean(recap_pred_Arb(:,3:4)),'--', 'LineWidth', 1.5, 'Color',[176 60 8]/255);
plotSpread([Behavior.GLME_OLunc.Reffects.PastTok_LowOLU Behavior.GLME_OLunc.Reffects.PastTok_HighOLU],...
    'xValues',[0.8 1.2], 'distributionColors',[79 184 255]/255,'spreadWidth',0.2);
plotSpread([Behavior.GLME_OLunc.Reffects.PastAct_LowOLU Behavior.GLME_OLunc.Reffects.PastAct_HighOLU],...
    'xValues',[1.8 2.2],'distributionColors',[246 146 100]/255,'spreadWidth',0.2);
errorbar([0.8 1.2],Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),Behavior.GLME_OLunc.Coefficients.SE([4 5]),...
    'Color',[21 160 255]/255,'LineWidth',1.5)
errorbar([1.8 2.2],Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),Behavior.GLME_OLunc.Coefficients.SE([2 3]),...
    'Color',[246 115 56]/255,'LineWidth',1.5)
errorbar([0.85 1.25],mean(recap_pred_Arb(:,1:2)),std(recap_pred_Arb(:,1:2))/sqrt(n_all),'.',...
    'Color',[0 114 189]/255,'LineWidth',1.5)
errorbar([1.85 2.25],mean(recap_pred_Arb(:,3:4)),std(recap_pred_Arb(:,3:4))/sqrt(n_all),'.',...
    'Color',[176 60 8]/255,'LineWidth',1.5)
xlim([0.50 2.4]); xticks([0.8 1.2 1.8 2.2]); xticklabels({'Low','High','Low','High'}); ylim([-1 2.5])
xlabel('OL uncertainty')
ylabel('ME-GLM effect')
legend({'Past outcome (EL effect) - data','Past action (OL effect) - data',...
     'Past outcome (EL effect) - model','Past action (OL effect) - model'})
set(gca,'box','off')

%GLM effect of EL uncertainty (trial definition)
recap_pred_EL = [PPC.GLME_ELunc.PastTok_LowELU(:,1) PPC.GLME_ELunc.PastTok_HighELU(:,1) ...
    PPC.GLME_ELunc.PastAct_LowELU(:,1) PPC.GLME_ELunc.PastAct_HighELU(:,1)];
recap_pred_OL = [PPC.GLME_ELunc.PastTok_LowELU(:,2) PPC.GLME_ELunc.PastTok_HighELU(:,2) ...
    PPC.GLME_ELunc.PastAct_LowELU(:,2) PPC.GLME_ELunc.PastAct_HighELU(:,2)];
recap_pred_Arb = [PPC.GLME_ELunc.PastTok_LowELU(:,3) PPC.GLME_ELunc.PastTok_HighELU(:,3) ...
    PPC.GLME_ELunc.PastAct_LowELU(:,3) PPC.GLME_ELunc.PastAct_HighELU(:,3)];
figure;
subplot(1,4,1); hold on
plot(0.85:1.85,Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5);
plot(1.15:2.15,Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_ELunc.Reffects.PastTok_LowELU Behavior.GLME_ELunc.Reffects.PastTok_HighELU],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_ELunc.Reffects.PastAct_LowELU Behavior.GLME_ELunc.Reffects.PastAct_HighELU],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),Behavior.GLME_ELunc.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),Behavior.GLME_ELunc.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1 2.5])
xlabel('EL uncertainty')
ylabel('ME-GLM effect')
title('Data')
subplot(1,4,2); hold on
plot(0.85:1.85,mean(recap_pred_EL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_EL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_EL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_EL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_EL(:,1:2)),std(recap_pred_EL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_EL(:,3:4)),std(recap_pred_EL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1 2.5])
xlabel('EL uncertainty')
title('EL model')
subplot(1,4,3); hold on
plot(0.85:1.85,mean(recap_pred_OL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_OL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_OL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_OL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_OL(:,1:2)),std(recap_pred_OL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_OL(:,3:4)),std(recap_pred_OL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1 2.5])
xlabel('EL uncertainty')
title('OL model')
subplot(1,4,4); hold on
plot(0.85:1.85,mean(recap_pred_Arb(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_Arb(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_Arb(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_Arb(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_Arb(:,1:2)),std(recap_pred_Arb(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_Arb(:,3:4)),std(recap_pred_Arb(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1 2.5])
xlabel('EL uncertainty')
title('DynArb model')
legend({'Past outcome (EL effect)','Past action (OL effect)'})

%Simplified plot (Figure 4F)
figure;
plot([0.8 1.2],Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5, 'Color',[21 160 255]/255); hold on
plot([1.8 2.2],Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5, 'Color',[246 115 56]/255);
plot([0.85 1.25],mean(recap_pred_Arb(:,1:2)),'--', 'LineWidth', 1.5, 'Color',[0 114 189]/255);
plot([1.85 2.25],mean(recap_pred_Arb(:,3:4)),'--', 'LineWidth', 1.5, 'Color',[176 60 8]/255);
plotSpread([Behavior.GLME_ELunc.Reffects.PastTok_LowELU Behavior.GLME_ELunc.Reffects.PastTok_HighELU],...
    'xValues',[0.8 1.2], 'distributionColors',[79 184 255]/255,'spreadWidth',0.2);
plotSpread([Behavior.GLME_ELunc.Reffects.PastAct_LowELU Behavior.GLME_ELunc.Reffects.PastAct_HighELU],...
    'xValues',[1.8 2.2],'distributionColors',[246 146 100]/255,'spreadWidth',0.2);
errorbar([0.8 1.2],Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),Behavior.GLME_ELunc.Coefficients.SE([4 5]),...
    'Color',[21 160 255]/255,'LineWidth',1.5)
errorbar([1.8 2.2],Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),Behavior.GLME_ELunc.Coefficients.SE([2 3]),...
    'Color',[246 115 56]/255,'LineWidth',1.5)
errorbar([0.85 1.25],mean(recap_pred_Arb(:,1:2)),std(recap_pred_Arb(:,1:2))/sqrt(n_all),'.',...
    'Color',[0 114 189]/255,'LineWidth',1.5)
errorbar([1.85 2.25],mean(recap_pred_Arb(:,3:4)),std(recap_pred_Arb(:,3:4))/sqrt(n_all),'.',...
    'Color',[176 60 8]/255,'LineWidth',1.5)
xlim([0.50 2.4]); xticks([0.8 1.2 1.8 2.2]); xticklabels({'Low','High','Low','High'}); ylim([-0.5 2.5])
xlabel('EL uncertainty')
ylabel('ME-GLM effect')
% legend({'Past outcome (EL effect) - data','Past action (OL effect) - data',...
%     'Past outcome (EL effect) - model','Past action (OL effect) - model'})
set(gca,'box','off')


%GLM effect of magnitude
recap_pred_EL = [PPC.GLME_Mag.PastTok_HighMag(:,1) PPC.GLME_Mag.PastTok_LowMag(:,1) ...
    PPC.GLME_Mag.PastAct_HighMag(:,1) PPC.GLME_Mag.PastAct_LowMag(:,1)];
recap_pred_OL = [PPC.GLME_Mag.PastTok_HighMag(:,2) PPC.GLME_Mag.PastTok_LowMag(:,2) ...
    PPC.GLME_Mag.PastAct_HighMag(:,2) PPC.GLME_Mag.PastAct_LowMag(:,2)];
recap_pred_Arb = [PPC.GLME_Mag.PastTok_HighMag(:,3) PPC.GLME_Mag.PastTok_LowMag(:,3) ...
    PPC.GLME_Mag.PastAct_HighMag(:,3) PPC.GLME_Mag.PastAct_LowMag(:,3)];
figure;
subplot(1,4,1); hold on
plot(0.85:1.85,Behavior.GLME_Mag.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5);
plot(1.15:2.15,Behavior.GLME_Mag.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_Mag.Reffects.PastTok_HighMag Behavior.GLME_Mag.Reffects.PastTok_LowMag],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_Mag.Reffects.PastAct_HighMag Behavior.GLME_Mag.Reffects.PastAct_LowMag],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_Mag.Coefficients.Estimate([4 5]),Behavior.GLME_Mag.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_Mag.Coefficients.Estimate([2 3]),Behavior.GLME_Mag.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'High','Low'}); ylim([-1 3])
xlabel('Magnitude')
ylabel('ME-GLM effect')
title('Data')
subplot(1,4,2); hold on
plot(0.85:1.85,mean(recap_pred_EL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_EL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_EL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_EL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_EL(:,1:2)),std(recap_pred_EL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_EL(:,3:4)),std(recap_pred_EL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'High','Low'}); ylim([-1 3])
xlabel('Magnitude')
title('EL model')
subplot(1,4,3); hold on
plot(0.85:1.85,mean(recap_pred_OL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_OL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_OL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_OL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_OL(:,1:2)),std(recap_pred_OL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_OL(:,3:4)),std(recap_pred_OL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'High','Low'}); ylim([-1 3])
xlabel('Magnitude')
title('OL model')
subplot(1,4,4); hold on
plot(0.85:1.85,mean(recap_pred_Arb(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_Arb(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_Arb(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_Arb(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_Arb(:,1:2)),std(recap_pred_Arb(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_Arb(:,3:4)),std(recap_pred_Arb(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'High','Low'}); ylim([-1 3])
xlabel('Magnitude')
title('DynArb model')
legend({'Past outcome (EL effect)','Past action (OL effect)'})


%% Correlation actual data and prediction of DynArb only (Figure S5A & S5C)
figure
subplot(2,3,1); hold on
diff_actual = Behavior.GLME_OLunc.Reffects.PastTok_LowOLU - Behavior.GLME_OLunc.Reffects.PastTok_HighOLU;
diff_sim = PPC.GLME_OLunc.PastTok_LowOLU(:,3) - PPC.GLME_OLunc.PastTok_HighOLU(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5);
lsline()
xlabel('actual'); ylabel('model'); xlim([-0.7 1]); ylim([-0.7 1])
title(["Low vs High OL unc",['R=',num2str(corr(diff_actual,diff_sim))]]);
subplot(2,3,2); hold on
diff_actual = Behavior.GLME_ELunc.Reffects.PastTok_LowELU - Behavior.GLME_ELunc.Reffects.PastTok_HighELU;
diff_sim = PPC.GLME_ELunc.PastTok_LowELU(:,3) - PPC.GLME_ELunc.PastTok_HighELU(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5);
lsline()
xlabel('actual'); ylabel('model'); xlim([-0.3 1.3]); ylim([-0.3 1.3])
title(["Low vs High EL unc",['R=',num2str(corr(diff_actual,diff_sim))]]);
subplot(2,3,3); hold on
diff_actual = Behavior.GLME_Mag.Reffects.PastTok_HighMag - Behavior.GLME_Mag.Reffects.PastTok_LowMag;
diff_sim = PPC.GLME_Mag.PastTok_HighMag(:,3) - PPC.GLME_Mag.PastTok_LowMag(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5);
lsline()
xlabel('actual'); ylabel('model'); xlim([-1 3]); ylim([-1 3])
title(["High vs Low Mag",['R=',num2str(corr(diff_actual,diff_sim))]]);
subplot(2,3,4); hold on
diff_actual = Behavior.GLME_OLunc.Reffects.PastAct_LowOLU - Behavior.GLME_OLunc.Reffects.PastAct_HighOLU;
diff_sim = PPC.GLME_OLunc.PastAct_LowOLU(:,3) - PPC.GLME_OLunc.PastAct_HighOLU(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor',OL_col);
lsline()
xlabel('actual'); ylabel('model'); xlim([-3.5 4]); ylim([-3.5 4])
title(['R=',num2str(corr(diff_actual,diff_sim))]);
subplot(2,3,5); hold on
diff_actual = Behavior.GLME_ELunc.Reffects.PastAct_LowELU - Behavior.GLME_ELunc.Reffects.PastAct_HighELU;
diff_sim = PPC.GLME_ELunc.PastAct_LowELU(:,3) - PPC.GLME_ELunc.PastAct_HighELU(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor',OL_col);
lsline()
xlabel('actual'); ylabel('model'); xlim([-0.6 0.5]); ylim([-0.6 0.5])
title(['R=',num2str(corr(diff_actual,diff_sim))]);
subplot(2,3,6); hold on
diff_actual = Behavior.GLME_Mag.Reffects.PastAct_HighMag - Behavior.GLME_Mag.Reffects.PastAct_LowMag;
diff_sim = PPC.GLME_Mag.PastAct_HighMag(:,3) - PPC.GLME_Mag.PastAct_LowMag(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor',OL_col);
lsline()
xlabel('actual'); ylabel('model'); xlim([-0.5 0.15]); ylim([-0.5 0.15])
title(['R=',num2str(corr(diff_actual,diff_sim))]);


%% repeat some of the above for block definitions 

%GLM effect of OL uncertainty (block definition)
recap_pred_EL = [PPC.GLME_OLunc_des.PastTok_LowOLU(:,1) PPC.GLME_OLunc_des.PastTok_HighOLU(:,1) ...
    PPC.GLME_OLunc_des.PastAct_LowOLU(:,1) PPC.GLME_OLunc_des.PastAct_HighOLU(:,1)];
recap_pred_OL = [PPC.GLME_OLunc_des.PastTok_LowOLU(:,2) PPC.GLME_OLunc_des.PastTok_HighOLU(:,2) ...
    PPC.GLME_OLunc_des.PastAct_LowOLU(:,2) PPC.GLME_OLunc_des.PastAct_HighOLU(:,2)];
recap_pred_Arb = [PPC.GLME_OLunc_des.PastTok_LowOLU(:,3) PPC.GLME_OLunc_des.PastTok_HighOLU(:,3) ...
    PPC.GLME_OLunc_des.PastAct_LowOLU(:,3) PPC.GLME_OLunc_des.PastAct_HighOLU(:,3)];
figure;
subplot(1,4,1); hold on
plot(0.85:1.85,Behavior.GLME_OLunc_des.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5);
plot(1.15:2.15,Behavior.GLME_OLunc_des.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_OLunc_des.Reffects.PastTok_LowOLU Behavior.GLME_OLunc_des.Reffects.PastTok_HighOLU],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_OLunc_des.Reffects.PastAct_LowOLU Behavior.GLME_OLunc_des.Reffects.PastAct_HighOLU],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_OLunc_des.Coefficients.Estimate([4 5]),Behavior.GLME_OLunc_des.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_OLunc_des.Coefficients.Estimate([2 3]),Behavior.GLME_OLunc_des.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1.5 3])
xlabel('OL uncertainty')
ylabel('ME-GLM effect')
title('Data')
subplot(1,4,2); hold on
plot(0.85:1.85,mean(recap_pred_EL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_EL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_EL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_EL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_EL(:,1:2)),std(recap_pred_EL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_EL(:,3:4)),std(recap_pred_EL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1.5 3])
xlabel('OL uncertainty')
title('EL model')
subplot(1,4,3); hold on
plot(0.85:1.85,mean(recap_pred_OL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_OL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_OL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_OL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_OL(:,1:2)),std(recap_pred_OL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_OL(:,3:4)),std(recap_pred_OL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1.5 3])
xlabel('OL uncertainty')
title('OL model')
subplot(1,4,4); hold on
plot(0.85:1.85,mean(recap_pred_Arb(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_Arb(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_Arb(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_Arb(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_Arb(:,1:2)),std(recap_pred_Arb(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_Arb(:,3:4)),std(recap_pred_Arb(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1.5 3])
xlabel('OL uncertainty')
title('DynArb model')
legend({'Past outcome (EL effect)','Past action (OL effect)'})

%GLM effect of EL uncertainty (block definition)
recap_pred_EL = [PPC.GLME_ELunc_des.PastTok_LowELU(:,1) PPC.GLME_ELunc_des.PastTok_HighELU(:,1) ...
    PPC.GLME_ELunc_des.PastAct_LowELU(:,1) PPC.GLME_ELunc_des.PastAct_HighELU(:,1)];
recap_pred_OL = [PPC.GLME_ELunc_des.PastTok_LowELU(:,2) PPC.GLME_ELunc_des.PastTok_HighELU(:,2) ...
    PPC.GLME_ELunc_des.PastAct_LowELU(:,2) PPC.GLME_ELunc_des.PastAct_HighELU(:,2)];
recap_pred_Arb = [PPC.GLME_ELunc_des.PastTok_LowELU(:,3) PPC.GLME_ELunc_des.PastTok_HighELU(:,3) ...
    PPC.GLME_ELunc_des.PastAct_LowELU(:,3) PPC.GLME_ELunc_des.PastAct_HighELU(:,3)];
figure;
subplot(1,4,1); hold on
plot(0.85:1.85,Behavior.GLME_ELunc_des.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5);
plot(1.15:2.15,Behavior.GLME_ELunc_des.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_ELunc_des.Reffects.PastTok_LowELU Behavior.GLME_ELunc_des.Reffects.PastTok_HighELU],...
    'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread([Behavior.GLME_ELunc_des.Reffects.PastAct_LowELU Behavior.GLME_ELunc_des.Reffects.PastAct_HighELU],...
    'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,Behavior.GLME_ELunc_des.Coefficients.Estimate([4 5]),Behavior.GLME_ELunc_des.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,Behavior.GLME_ELunc_des.Coefficients.Estimate([2 3]),Behavior.GLME_ELunc_des.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1 2.5])
xlabel('EL uncertainty')
ylabel('ME-GLM effect')
title('Data')
subplot(1,4,2); hold on
plot(0.85:1.85,mean(recap_pred_EL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_EL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_EL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_EL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_EL(:,1:2)),std(recap_pred_EL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_EL(:,3:4)),std(recap_pred_EL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1 2.5])
xlabel('EL uncertainty')
title('EL model')
subplot(1,4,3); hold on
plot(0.85:1.85,mean(recap_pred_OL(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_OL(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_OL(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_OL(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_OL(:,1:2)),std(recap_pred_OL(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_OL(:,3:4)),std(recap_pred_OL(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1 2.5])
xlabel('EL uncertainty')
title('OL model')
subplot(1,4,4); hold on
plot(0.85:1.85,mean(recap_pred_Arb(:,1:2)),'-', 'LineWidth', 1.5);
plot(1.15:2.15,mean(recap_pred_Arb(:,3:4)),'-', 'LineWidth', 1.5);
plotSpread(recap_pred_Arb(:,1:2),'xValues',(1:2)-0.15,'distributionColors',[122/256 165/256 207/256],'spreadWidth',0.2);
plotSpread(recap_pred_Arb(:,3:4),'xValues',(1:2)+0.15,'distributionColors',[214/256 139/256 135/256],'spreadWidth',0.2);
errorbar(0.85:1.85,mean(recap_pred_Arb(:,1:2)),std(recap_pred_Arb(:,1:2))/sqrt(n_all),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1.15:2.15,mean(recap_pred_Arb(:,3:4)),std(recap_pred_Arb(:,3:4))/sqrt(n_all),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35]); xticks([1 2]); xticklabels({'Low','High'}); ylim([-1 2.5])
xlabel('EL uncertainty')
title('DynArb model')
legend({'Past outcome (EL effect)','Past action (OL effect)'})


%% Correlation actual data and prediction of DynArb only (Figure S5E & S5G)
figure
subplot(2,2,1); hold on
diff_actual = Behavior.GLME_OLunc_des.Reffects.PastTok_LowOLU - Behavior.GLME_OLunc_des.Reffects.PastTok_HighOLU;
diff_sim = PPC.GLME_OLunc_des.PastTok_LowOLU(:,3) - PPC.GLME_OLunc_des.PastTok_HighOLU(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5);
lsline()
xlabel('actual'); ylabel('model'); xlim([-0.5 0.5]); ylim([-0.5 0.5])
title(["Low vs High OL unc",['R=',num2str(corr(diff_actual,diff_sim))]]);
subplot(2,2,2); hold on
diff_actual = Behavior.GLME_ELunc_des.Reffects.PastTok_LowELU - Behavior.GLME_ELunc_des.Reffects.PastTok_HighELU;
diff_sim = PPC.GLME_ELunc_des.PastTok_LowELU(:,3) - PPC.GLME_ELunc_des.PastTok_HighELU(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5);
lsline()
xlabel('actual'); ylabel('model'); xlim([-0.3 0.7]); ylim([-0.3 0.7])
title(["Low vs High EL unc",['R=',num2str(corr(diff_actual,diff_sim))]]);
subplot(2,2,3); hold on
diff_actual = Behavior.GLME_OLunc_des.Reffects.PastAct_LowOLU - Behavior.GLME_OLunc_des.Reffects.PastAct_HighOLU;
diff_sim = PPC.GLME_OLunc_des.PastAct_LowOLU(:,3) - PPC.GLME_OLunc_des.PastAct_HighOLU(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor',OL_col);
lsline()
xlabel('actual'); ylabel('model'); xlim([-0.5 1]); ylim([-0.5 1])
title(['R=',num2str(corr(diff_actual,diff_sim))]);
subplot(2,2,4); hold on
diff_actual = Behavior.GLME_ELunc_des.Reffects.PastAct_LowELU - Behavior.GLME_ELunc_des.Reffects.PastAct_HighELU;
diff_sim = PPC.GLME_ELunc_des.PastAct_LowELU(:,3) - PPC.GLME_ELunc_des.PastAct_HighELU(:,3);
scatter(diff_actual,diff_sim,'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor',OL_col);
lsline()
xlabel('actual'); ylabel('model'); xlim([-0.6 0.4]); ylim([-0.6 0.4])
title(['R=',num2str(corr(diff_actual,diff_sim))]);
