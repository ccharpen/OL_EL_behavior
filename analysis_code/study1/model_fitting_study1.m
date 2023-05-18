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
data_all = cell(n_all,1);
for s=1:n_all
    subNb = subID_list(s);
    data_all{s} = table2array(data(data.subNb==subNb,2:end));
end

%specify loglikelihood functions
func_list = {@LL_Baseline; @LL_ExpLearn; @LL_ObsLearn; @LL_FixArb; @LL_DynArb};
n_mod = length(func_list);

%specify output files
out_dir = 'model_fitting_outputs';
if ~exist(out_dir, 'dir')
   mkdir(out_dir)
end
out_fname_list = {'lap_Baseline.mat';'lap_ExpLearn.mat';'lap_ObsLearn.mat';...
    'lap_FixArb.mat';'lap_DynArb.mat'};

mod_names = {'Baseline'; 'ExpLearn'; 'ObsLearn'; 'FixArb'; 'DynArb'};

%number of parameters for each model
np = [4;3;2;6;6]; 

%specify parameter priors
v = 6.25; %parameter variance (6.25 is large enough to cover a wide range of parameters with no excessive penalty)

%% Individual-level fits
%run cbm_lap for each model
%this will fit every model to each subject's data separately (ie in a
%non-hierarchical fashion), using Laplace approximation, which needs a
%normal prior for every parameter
numcores = feature('numcores');
parfor (m=1:n_mod,numcores-1) 
    prior = struct('mean',zeros(np(m),1),'variance',v);
    cbm_lap(data_all, func_list{m}, prior, [out_dir fs out_fname_list{m}]);
end

%transform parameters and calculate model fitting metrics
fitRecap = struct();
fitRecap.pseudoR2 = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
fitRecap.AIC = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
fitRecap.BIC = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
fitRecap.LHsub = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
fitRecap.corrsub = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
for m=1:n_mod
    fname = [out_dir fs out_fname_list{m}];
    load(fname,'cbm')
    params = cbm.output.parameters;
    npar = np(m);
    
    cbm.output.paramTrans = params;
    if m>=2
        cbm.output.paramTrans(:,1) = exp(params(:,1));
    end
    if m==2 || m==3 %single strategy model (one learning rate)
        cbm.output.paramTrans(:,2) = 1./(1+exp(-params(:,2)));        
    elseif m>=4 %arbitration models (second beta + 2 learning rates)    
        cbm.output.paramTrans(:,2) = exp(params(:,2));
        cbm.output.paramTrans(:,3) = 1./(1+exp(-params(:,3)));
        cbm.output.paramTrans(:,4) = 1./(1+exp(-params(:,4)));
    end
    if m==4 %FixArb model (weight from 0-1)
        cbm.output.paramTrans(:,5) = 1./(1+exp(-params(:,5)));
    end
    
    for s=1:n_all
        subNb = subID_list(s);
        P = table2array(data(data.subNb==subNb,2:end));
        ntg = sum(P(:,19)==0);
        [ll,P_pred] = func_list{m}(params(s,:),P);
        %calculate model
        cbm.output.pseudoR2(s,1) = 1 - ll/(ntg*log(0.5));
        cbm.output.BIC(s,1) = -2*ll + log(ntg)*npar;
        cbm.output.AIC(s,1) = -2*ll + 2*npar;
        if cbm.output.pseudoR2(s,1)<0
            cbm.output.pseudoR2(s,1)=0;
        end
        %using values generated by LL function
        cbm.output.LHsub(s,1) = nanmean(P_pred(:,1));
        cbm.output.corrsub(s,1) = nanmean(P_pred(:,3));   
    end
    save(fname,'cbm')
    
    fitRecap.pseudoR2{:,m} = cbm.output.pseudoR2;
    fitRecap.AIC{:,m}      = cbm.output.AIC;
    fitRecap.BIC{:,m}      = cbm.output.BIC;
    fitRecap.LHsub{:,m}    = cbm.output.LHsub;
    fitRecap.corrsub{:,m}  = cbm.output.corrsub;
    
    fitRecap.paramRaw.(mod_names{m}) = cbm.output.parameters;
    fitRecap.paramTrans.(mod_names{m}) = cbm.output.paramTrans;
    
end
save([out_dir fs 'Recap_model_fitting.mat'],'fitRecap');

%% Hierarchical fits on single models
fname_hbi = {'hbi_Baseline.mat';'hbi_ExpLearn.mat';'hbi_ObsLearn.mat';'hbi_FixArb.mat';'hbi_DynArb.mat'};
parfor (m=1:n_mod,numcores-1) 
    cbm_hbi(data_all, func_list(m), [out_dir fs out_fname_list(m)], [out_dir fs fname_hbi{m}]);
end

%calculate AIC and other measures based on hierarchical fits
hfitRecap = struct();
hfitRecap.pseudoR2 = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
hfitRecap.AIC = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
hfitRecap.BIC = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
hfitRecap.LHsub = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
hfitRecap.corrsub = table('Size',[n_all n_mod],'VariableTypes',repmat({'double'},1,n_mod),'VariableNames',mod_names);
for m=1:n_mod
    fname = [out_dir fs fname_hbi{m}];
    load(fname,'cbm')
    params = cbm.output.parameters{:};
    npar = np(m);
    
    cbm.output.paramTrans = params;
    if m>=2
        cbm.output.paramTrans(:,1) = exp(params(:,1));
    end
    if m==2 || m==3 %single strategy model (one learning rate)
        cbm.output.paramTrans(:,2) = 1./(1+exp(-params(:,2)));        
    elseif m>=4 %arbitration models (second beta + 2 learning rates)    
        cbm.output.paramTrans(:,2) = exp(params(:,2));
        cbm.output.paramTrans(:,3) = 1./(1+exp(-params(:,3)));
        cbm.output.paramTrans(:,4) = 1./(1+exp(-params(:,4)));
    end
    if m==4 %FixArb model (weight from 0-1)
        cbm.output.paramTrans(:,5) = 1./(1+exp(-params(:,5)));
    end
    
    for s=1:n_all
        subNb = subID_list(s);
        P = table2array(data(data.subNb==subNb,2:end));
        ntg = sum(P(:,19)==0);
        [ll,P_pred] = func_list{m}(params(s,:),P);
        %calculate model
        cbm.output.pseudoR2(s,1) = 1 - ll/(ntg*log(0.5));
        cbm.output.BIC(s,1) = -2*ll + log(ntg)*npar;
        cbm.output.AIC(s,1) = -2*ll + 2*npar;
        if cbm.output.pseudoR2(s,1)<0
            cbm.output.pseudoR2(s,1)=0;
        end
        %using values generated by LL function
        cbm.output.LHsub(s,1) = nanmean(P_pred(:,1));
        cbm.output.corrsub(s,1) = nanmean(P_pred(:,3));   
    end
    save(fname,'cbm')
    
    hfitRecap.pseudoR2{:,m} = cbm.output.pseudoR2;
    hfitRecap.AIC{:,m}      = cbm.output.AIC;
    hfitRecap.BIC{:,m}      = cbm.output.BIC;
    hfitRecap.LHsub{:,m}    = cbm.output.LHsub;
    hfitRecap.corrsub{:,m}  = cbm.output.corrsub;
    
    hfitRecap.paramRaw.(mod_names{m}) = cbm.output.parameters;
    hfitRecap.paramTrans.(mod_names{m}) = cbm.output.paramTrans;
end
save([out_dir fs 'Recap_model_fitting.mat'],'fitRecap','hfitRecap');

%% parameter recovery
%generate data for each model using each subject's parameters
n_sim = 10;
gc_mod_list = {@generate_choice_Baseline; @generate_choice_ExpLearn; ...
    @generate_choice_ObsLearn; @generate_choice_FixArb; @generate_choice_DynArb};
load([out_dir fs 'Recap_model_fitting.mat'])
Actual_params = fitRecap.paramRaw;
%create result directory if it doesn't exist
if ~exist('parameter_recovery', 'dir')
   mkdir('parameter_recovery')
end
for m = 1:n_mod
    m
    params_mod = Actual_params.(mod_names{m});
    npar = np(m);
    parfor (n=1:n_sim,5)
        data_sim_all = cell(n_all,1);
        for s = 1:n_all
            subNb = subID_list(s);
            P = table2array(data(data.subNb==subNb,2:end));
            param_sub = params_mod(s,:);
            P_pred = gc_mod_list{m}(param_sub,P);
            P_new = P;
            P_new(:,14) = P_pred(:,4); %choice (1: orange, 0: blue)
            P_new(:,17) = P_pred(:,5); %token  (1: orange, 2: blue)
            P_new(:,16) = P_pred(:,6); %whether predicted choice is correct (1) or not (0)
            P_new(:,18) = P_pred(:,8)*100; %predicted outcome, should be scaled from 0 to 100
            if m==1
                P_new(:,13) = P_pred(:,9); %action (1: left, 0: right)
                %not needed for other models than baseline and integrated
            end
            data_sim_all{s,1} = P_new;
        end
        prior = struct('mean',zeros(npar,1),'variance',v);
        out_fname = ['parameter_recovery' fs 'lap_' mod_names{m} '_sim' num2str(n) '.mat'];
        cbm_lap(data_sim_all, func_list{m}, prior, out_fname);
        out_fname_hbi = ['parameter_recovery' fs 'hbi_' mod_names{m} '_sim' num2str(n) '.mat'];
        cbm_hbi(data_sim_all, func_list(m), {out_fname}, {out_fname_hbi});
    end
end

%now load the cbm files and compute correlations between actual and
%recovered parameters
for m=1:n_mod
    npar = np(m);
    recov_params_indiv = zeros(n_all,npar,n_sim);
    recov_params_hbi = zeros(n_all,npar,n_sim);
    for n=1:n_sim
        out_fname = ['parameter_recovery' fs 'lap_' mod_names{m} '_sim' num2str(n) '.mat'];
        res = load(out_fname);
        recov_params_indiv(:,:,n) = res.cbm.output.parameters;
        out_fname_hbi = ['parameter_recovery' fs 'hbi_' mod_names{m} '_sim' num2str(n) '.mat'];
        res = load(out_fname_hbi);
        recov_params_hbi(:,:,n) = cell2mat(res.cbm.output.parameters);
    end
    Recov_params_indiv.(mod_names{m}) = mean(recov_params_indiv,3);
    Recov_params_hbi.(mod_names{m}) = mean(recov_params_indiv,3);
    for p=1:npar
        CorrelActRecov_indiv.(mod_names{m})(1,p) = corr(Actual_params.(mod_names{m})(:,p),Recov_params_indiv.(mod_names{m})(:,p));
        CorrelActRecov_hbi.(mod_names{m})(1,p) = corr(Actual_params.(mod_names{m})(:,p),Recov_params_hbi.(mod_names{m})(:,p));
    end
end
save(['parameter_recovery' fs 'Recap_parameter_recovery.mat'],'Actual_params','Recov_params_indiv',...
    'Recov_params_hbi','CorrelActRecov_indiv','CorrelActRecov_hbi');

%% plots for parameter recovery
load(['parameter_recovery' fs 'Recap_parameter_recovery.mat'],'Actual_params','Recov_params_hbi','CorrelActRecov_hbi');

p_cell = {{'ColorBias','HandBias','StickyAct','ActImit'}; {'EL beta','EL alpha','MagBoost'}; {'OL beta','OL alpha'}; ...
    {'OL beta','EL beta','OL alpha','EL alpha','w(OL>EL)','MagBoost'}; ...
    {'OL beta','EL beta','OL alpha','EL alpha','MagBoost','bias(OL>EL)'}};
for m=1:n_mod
    npar = np(m);
    par_rec_mat = zeros(npar,npar);
    ap = Actual_params.(mod_names{m});
    rp = Recov_params_hbi.(mod_names{m});
    p_list = p_cell{m};
    for pa = 1:npar
        for pr = 1:npar
            par_rec_mat(pa,pr) = corr(ap(:,pa),rp(:,pr));
        end
    end
    figure;
    h = heatmap(p_list,p_list,par_rec_mat);
    h.Title = ['Parameter Recovery (' mod_names{m} ' model)'];
    h.XLabel = 'Actual parameters';
    h.YLabel = 'Recovered parameters';
    h.CellLabelFormat = '%.3f';
end

%baseline model
act_par = Actual_params.Baseline;
rec_par = Recov_params_hbi.Baseline;
figure;
subplot(2,2,1); hold on
plot(act_par(:,1),rec_par(:,1),'.k');lsline
plot([min(act_par(:,1)) max(act_par(:,1))],[min(act_par(:,1)) max(act_par(:,1))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Color Bias')
subplot(2,2,2); hold on
plot(act_par(:,2),rec_par(:,2),'.k');lsline
plot([min(act_par(:,2)) max(act_par(:,2))],[min(act_par(:,2)) max(act_par(:,2))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Hand Bias')
subplot(2,2,3); hold on
plot(act_par(:,3),rec_par(:,3),'.k');lsline
plot([min(act_par(:,3)) max(act_par(:,3))],[min(act_par(:,3)) max(act_par(:,3))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Sticky Action')
subplot(2,2,4); hold on
plot(act_par(:,4),rec_par(:,4),'.k');lsline
plot([min(act_par(:,4)) max(act_par(:,4))],[min(act_par(:,4)) max(act_par(:,4))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Action Imitation')

%EL model
act_par = Actual_params.ExpLearn;
rec_par = Recov_params_hbi.ExpLearn;
figure;
subplot(1,3,1); hold on
plot(act_par(:,1),rec_par(:,1),'.k');lsline
plot([min(act_par(:,1)) max(act_par(:,1))],[min(act_par(:,1)) max(act_par(:,1))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Softmax beta')
subplot(1,3,2); hold on
plot(act_par(:,2),rec_par(:,2),'.k');lsline
plot([min(act_par(:,2)) max(act_par(:,2))],[min(act_par(:,2)) max(act_par(:,2))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('EL learning rate')
subplot(1,3,3); hold on
plot(act_par(:,3),rec_par(:,3),'.k');lsline
plot([min(act_par(:,3)) max(act_par(:,3))],[min(act_par(:,3)) max(act_par(:,3))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Magnitude boosting')

%OL model
act_par = Actual_params.ObsLearn;
rec_par = Recov_params_hbi.ObsLearn;
figure;
subplot(1,2,1); hold on
plot(act_par(:,1),rec_par(:,1),'.k');lsline;hold on
plot([min(act_par(:,1)) max(act_par(:,1))],[min(act_par(:,1)) max(act_par(:,1))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Softmax beta')
subplot(1,2,2); hold on
plot(act_par(:,2),rec_par(:,2),'.k');lsline;hold on
plot([min(act_par(:,2)) max(act_par(:,2))],[min(act_par(:,2)) max(act_par(:,2))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('OL learning rate')

%FixArb model
act_par = Actual_params.FixArb;
rec_par = Recov_params_hbi.FixArb;
figure;
subplot(2,3,1); hold on
plot(act_par(:,1),rec_par(:,1),'.k');lsline
plot([min(act_par(:,1)) max(act_par(:,1))],[min(act_par(:,1)) max(act_par(:,1))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('OL softmax beta')
subplot(2,3,2); hold on
plot(act_par(:,2),rec_par(:,2),'.k');lsline
plot([min(act_par(:,2)) max(act_par(:,2))],[min(act_par(:,2)) max(act_par(:,2))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('EL softmax beta')
subplot(2,3,3); hold on
plot(act_par(:,3),rec_par(:,3),'.k');lsline
plot([min(act_par(:,3)) max(act_par(:,3))],[min(act_par(:,3)) max(act_par(:,3))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('OL learning rate')
subplot(2,3,4); hold on
plot(act_par(:,4),rec_par(:,4),'.k');lsline
plot([min(act_par(:,4)) max(act_par(:,4))],[min(act_par(:,4)) max(act_par(:,4))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('EL learning rate')
subplot(2,3,5); hold on
plot(act_par(:,5),rec_par(:,5),'.k');lsline
plot([min(act_par(:,5)) max(act_par(:,5))],[min(act_par(:,5)) max(act_par(:,5))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Weight(OL>EL)')
subplot(2,3,6); hold on
plot(act_par(:,6),rec_par(:,6),'.k');lsline
plot([min(act_par(:,6)) max(act_par(:,6))],[min(act_par(:,6)) max(act_par(:,6))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Magnitude boosting')

%DynArb model
act_par = Actual_params.DynArb;
rec_par = Recov_params_hbi.DynArb;
figure;
subplot(2,3,1); hold on
plot(act_par(:,1),rec_par(:,1),'.k');lsline
plot([min(act_par(:,1)) max(act_par(:,1))],[min(act_par(:,1)) max(act_par(:,1))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('OL softmax beta')
subplot(2,3,2); hold on
plot(act_par(:,2),rec_par(:,2),'.k');lsline
plot([min(act_par(:,2)) max(act_par(:,2))],[min(act_par(:,2)) max(act_par(:,2))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('EL softmax beta')
subplot(2,3,3); hold on
plot(act_par(:,3),rec_par(:,3),'.k');lsline
plot([min(act_par(:,3)) max(act_par(:,3))],[min(act_par(:,3)) max(act_par(:,3))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('OL learning rate')
subplot(2,3,4); hold on
plot(act_par(:,4),rec_par(:,4),'.k');lsline
plot([min(act_par(:,4)) max(act_par(:,4))],[min(act_par(:,4)) max(act_par(:,4))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('EL learning rate')
subplot(2,3,5); hold on
plot(act_par(:,5),rec_par(:,5),'.k');lsline
plot([min(act_par(:,5)) max(act_par(:,5))],[min(act_par(:,5)) max(act_par(:,5))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Magnitude boosting')
subplot(2,3,6); hold on
plot(act_par(:,6),rec_par(:,6),'.k');lsline
plot([min(act_par(:,6)) max(act_par(:,6))],[min(act_par(:,6)) max(act_par(:,6))],'--r');
xlabel('Actual'); ylabel('Recovered'); title('Bias(OL>EL)')


%% Hierarchical fitting across all 5 models
models = {@LL_Baseline; @LL_ExpLearn; @LL_ObsLearn; @LL_FixArb; @LL_DynArb};
fcbm_maps = {'lap_Baseline.mat';'lap_ExpLearn.mat';'lap_ObsLearn.mat';...
    'lap_FixArb.mat';'lap_DynArb.mat'};
fname_hbi = 'hbi_5mods.mat';
cd(out_dir)
cbm_hbi(data_all, models, fcbm_maps, fname_hbi);
cd ..

%% check that OL vs EL models predict corresponding behavioral signature
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


%% plot trial-by-trial arbitration weight for an example subject with bias~0
i_small_bias = abs(fitRecap.paramRaw.DynArb(:,6))<0.1 & fitRecap.paramRaw.DynArb(:,1)>1 & fitRecap.paramRaw.DynArb(:,2)>1;
params_dynarb = fitRecap.paramRaw.DynArb(i_small_bias,:);
subNb = subID_list(find(i_small_bias));
P = table2array(data(data.subNb==subNb,2:end));
P_pred = generate_choice_DynArb(params_dynarb, P);
r_OL = P_pred(:,23);
r_EL = P_pred(:,24);
w = P_pred(:,25);
OL_lu = zeros(ntr,1);
EL_lu = zeros(ntr,1);
RM_h = zeros(ntr,1);
outc_bin = P(:,18); outc_bin(outc_bin>0)=1;
for t=1:ntr
    if P(t,4)>=3
        if P(t-1,6)==P(t,6)
            OL_lu(t) = 1;
        end
        if (outc_bin(t-2)==outc_bin(t-1) && P(t-2,14)==P(t-1,14)) || (outc_bin(t-2)~=outc_bin(t-1) && P(t-2,14)~=P(t-1,14))
            %same choice, same outcome OR different choice, different outcome
            EL_lu(t) = 1;
        end
    end
    if P(t,3)>=2 && P(t-1,18)>25
        RM_h(t) = 1;
    end
end

figure;
plot(~EL_lu*0.85+0.09,'.','Color','#0072BD','MarkerSize',8); hold on
plot(OL_lu*0.91+0.06,'.','Color','#D95319','MarkerSize',8);
plot(~RM_h*0.97+0.03,'.','Color','#EDB120','MarkerSize',8);
plot(w,'k','LineWidth',2)
ylim([0 1])
set(gca,'box','off')
xlabel('Trial')
legend({'EL uncertainty','OL uncertainty','Reward Magnitude','Arbitration weight'})

%% extract trial-by-trial arbitration weight and plot mean per condition
ntr = 160;
w_bd = nan(n_all, 8);
for s=1:n_all
    params_dynarb = fitRecap.paramRaw.DynArb(s,:);
    subNb = subID_list(s);
    P = table2array(data(data.subNb==subNb,2:end));
    P_pred = generate_choice_DynArb(params_dynarb, P);
    w = P_pred(:,25);
    
    %define uncertainty condition based on key trials, instead of design conditions
    %OL unc = low if past 2 box-token transition consistent; high if not
    %EL unc = low if W-st-W, W-sw-0, 0-st-0, 0-sw-W; high otherwise
    %RM = high if >25, low otherwise
    OL_lu = zeros(ntr,1);
    EL_lu = zeros(ntr,1);
    RM_h = zeros(ntr,1);
    outc_bin = P(:,18); outc_bin(outc_bin>0)=1;
    for t=1:ntr
        if P(t,4)>=3
            if P(t-1,6)==P(t,6)
                OL_lu(t) = 1;
            end
            if (outc_bin(t-2)==outc_bin(t-1) && P(t-2,14)==P(t-1,14)) || (outc_bin(t-2)~=outc_bin(t-1) && P(t-2,14)~=P(t-1,14))
                %same choice, same outcome OR different choice, different outcome
                EL_lu(t) = 1;
            end
        end
        if P(t,3)>=2 && P(t-1,18)>25
            RM_h(t) = 1;
        end
    end
    OL_lu = logical(OL_lu);
    EL_lu = logical(EL_lu);
    RM_h = logical(RM_h);
    
    w_bd(s,:) = [nanmean(w(OL_lu & EL_lu & RM_h)) nanmean(w(OL_lu & EL_lu & ~RM_h)) ...
        nanmean(w(OL_lu & ~EL_lu & RM_h)) nanmean(w(OL_lu & ~EL_lu & ~RM_h)) ...
        nanmean(w(~OL_lu & EL_lu & RM_h)) nanmean(w(~OL_lu & EL_lu & ~RM_h)) ...
        nanmean(w(~OL_lu & ~EL_lu & RM_h)) nanmean(w(~OL_lu & ~EL_lu & ~RM_h))];  
end

figure;
subplot(1,2,1); hold on
b = bar([nanmean(w_bd(:,1:2)); nanmean(w_bd(:,3:4))],0.8,'FaceColor','flat','EdgeColor','k','LineWidth',1);
b(1).CData = [0.9290 0.6940 0.1250]; 
b(2).CData = [0.4940 0.1840 0.5560];
plotSpread(w_bd(:,1:4),'xValues',[0.85 1.15 1.85 2.15],'distributionColors',[0.4 0.4 0.4]); 
errorbar([0.85 1.15 1.85 2.15],nanmean(w_bd(:,1:4)),nanstd(w_bd(:,1:4))/sqrt(n_all),'.k','LineWidth',1.5)
plot([0.5 2.5],[0.5 0.5],'--k')
xticks([1 2])
xticklabels({'Low','High'})
xlabel('EL uncertainty')
ylabel('Arbitration weight (OL vs EL)')
title('LOW OL uncertainty')
subplot(1,2,2); hold on
b = bar([nanmean(w_bd(:,5:6)); nanmean(w_bd(:,7:8))],0.8,'FaceColor','flat','EdgeColor','k','LineWidth',1);
b(1).CData = [0.9290 0.6940 0.1250]; 
b(2).CData = [0.4940 0.1840 0.5560];
plotSpread(w_bd(:,5:8),'xValues',[0.85 1.15 1.85 2.15],'distributionColors',[0.4 0.4 0.4]); 
errorbar([0.85 1.15 1.85 2.15],nanmean(w_bd(:,5:8)),nanstd(w_bd(:,5:8))/sqrt(n_all),'.k','LineWidth',1.5)
plot([0.5 2.5],[0.5 0.5],'--k')
xticks([1 2])
xticklabels({'Low','High'})
xlabel('EL uncertainty')
title('HIGH OL uncertainty')
leg = legend({'High','Low'});
title(leg,'Magnitude')
