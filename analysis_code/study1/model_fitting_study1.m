%This script performs all the model-fitting analyses for Study 1, parameter
%recovery analysis, and generates the results presented in Table 1 (AIC, Frequency, Nbest),
%Table S3 (AIC), as well as the plots presented in: Figures S3B-F and S6A-B

clear all
close all
fs = filesep;

data_dir = ['..' fs '..' fs 'data'];
mod_dir = ['..' fs 'dependencies' fs 'model_functions']; %directory where modelling functions are saved (common to both studies)

addpath(['..' fs 'dependencies']);
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

%AIC values reported in Table 1:
AIC_values = mean(table2array(fitRecap.AIC))

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

%now load the cbm files and compute correlations between actual and recovered parameters
for m=1:n_mod
    npar = np(m);
    actual_params = Actual_params.(mod_names{m});
    recap_corr_indiv = zeros(npar,npar,n_sim);
    recap_corr_hbi = zeros(npar,npar,n_sim);
    for n=1:n_sim
        out_fname = ['parameter_recovery' fs 'lap_' mod_names{m} '_sim' num2str(n) '.mat'];
        res = load(out_fname);
        params_indiv_sim = res.cbm.output.parameters;
        out_fname_hbi = ['parameter_recovery' fs 'hbi_' mod_names{m} '_sim' num2str(n) '.mat'];
        res = load(out_fname_hbi);
        params_hbi_sim = cell2mat(res.cbm.output.parameters);
        for pa=1:npar
            for pr=1:npar
                recap_corr_indiv(pa,pr,n) = corr(actual_params(:,pa),params_indiv_sim(:,pr));
                recap_corr_hbi(pa,pr,n) = corr(actual_params(:,pa),params_hbi_sim(:,pr));
            end
        end
    end
    Recovery_indiv.(mod_names{m}) = mean(recap_corr_indiv,3);
    Recovery_hbi.(mod_names{m}) = mean(recap_corr_hbi,3);
end
save(['parameter_recovery' fs 'Recap_parameter_recovery_revised.mat'],'Actual_params','Recovery_indiv','Recovery_hbi');

%% Plots for parameter recovery (Figures S3B-F)
load(['parameter_recovery' fs 'Recap_parameter_recovery_revised.mat'],'Actual_params','Recovery_indiv','Recovery_hbi');

p_cell = {{'ColorBias','HandBias','StickyAct','ActImit'}; {'EL beta','EL alpha','MagBoost'}; {'OL beta','OL alpha'}; ...
    {'OL beta','EL beta','OL alpha','EL alpha','w(OL>EL)','MagBoost'}; ...
    {'OL beta','EL beta','OL alpha','EL alpha','MagBoost','bias(OL>EL)'}};
for m=1:n_mod
    p_list = p_cell{m};
    rp = Recovery_indiv.(mod_names{m});
    figure;
    h = heatmap(p_list,p_list,rp);
    h.Title = ['Parameter Recovery (' mod_names{m} ' model)'];
    h.XLabel = 'Actual parameters';
    h.YLabel = 'Recovered parameters';
    h.CellLabelFormat = '%.3f';
    h.Colormap = redwhiteblue(-1,1);
    h.ColorLimits = [-1 1];
end

%% Hierarchical fitting across all 5 models
models = {@LL_Baseline; @LL_ExpLearn; @LL_ObsLearn; @LL_FixArb; @LL_DynArb};
fcbm_maps = {'lap_Baseline.mat';'lap_ExpLearn.mat';'lap_ObsLearn.mat';...
    'lap_FixArb.mat';'lap_DynArb.mat'};
fname_hbi = 'hbi_5mods.mat';
cd(out_dir)
cbm_hbi(data_all, models, fcbm_maps, fname_hbi);
cd ..
%the model frequency values reported in Table 1 can be found in
%hbi_5mods.mat output file, in variable cbm.output.model_frequency

%define groups based on hierarchical fit model responsibility values
load(['model_fitting_outputs' fs 'hbi_5mods.mat'])
hfit_recap = cbm.output.responsibility;
for s=1:n_all
    hfit_recap(s,6) = find(hfit_recap(s,1:5)==max(hfit_recap(s,1:5)));
end
group = hfit_recap(:,6);
gsize = [sum(group==1) sum(group==2) sum(group==3) sum(group==4) sum(group==5)];
gdist = gsize/sum(gsize);
save('Recap_model_fitting.mat','fitRecap','hfitRecap','group');
%gsize and gdist contain the values reported in the Nbest column of Table 1

%% Plot trial-by-trial arbitration weight for an example subject with bias~0 (Figure S6A)
i_small_bias = abs(fitRecap.paramRaw.DynArb(:,6))<0.1 & fitRecap.paramRaw.DynArb(:,1)>1 & fitRecap.paramRaw.DynArb(:,2)>1;
params_dynarb = fitRecap.paramRaw.DynArb(i_small_bias,:);
subNb = subID_list(find(i_small_bias));
P = table2array(data(data.subNb==subNb,2:end));
P_pred = generate_choice_DynArb(params_dynarb, P);
r_OL = P_pred(:,23);
r_EL = P_pred(:,24);
w = P_pred(:,25);
ntr=160;
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

%% Extract trial-by-trial arbitration weight and plot mean per condition (Figure S6B)
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


%% Model-fitting for additional EL mechanisms/learning of magnitude

%specify loglikelihood functions
func_list_EL = {@LL_ExpLearn;@LL_ExpLearn_nomag; @LL_ExpLearn_mag; @LL_ExpLearn_decay};
n_mod_EL = length(func_list_EL);

%specify output files
out_fname_list_EL = {'lap_ExpLearn.mat';'lap_ExpLearn_nomag.mat';'lap_ExpLearn_mag.mat';'lap_ExpLearn_decay.mat'};
mod_names_EL = {'ExpLearn';'ExpLearn_nomag'; 'ExpLearn_mag'; 'ExpLearn_decay'};

%number of parameters for each model
np = [3;2;2;4]; 

%specify parameter priors
v = 6.25; %parameter variance (6.25 is large enough to cover a wide range of parameters with no excessive penalty)

%individual model-fitting
parfor (m=2:n_mod_EL,3) %starts at 2 because no need to refit ExpLearn
    prior = struct('mean',zeros(np(m),1),'variance',v);
    cbm_lap(data_all, func_list_EL{m}, prior, [out_dir fs out_fname_list_EL{m}]);
end

%transform parameters and calculate model fitting metrics
fitRecapEL = struct();
fitRecapEL.pseudoR2 = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
fitRecapEL.AIC = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
fitRecapEL.BIC = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
fitRecapEL.LHsub = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
fitRecapEL.corrsub = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
for m=1:n_mod_EL
    fname = [out_dir fs out_fname_list_EL{m}];
    load(fname,'cbm')
    params = cbm.output.parameters;
    npar = np(m);
    
    if m >= 2
        cbm.output.paramTrans = params;
        cbm.output.paramTrans(:,1) = exp(params(:,1)); %softmax beta
        cbm.output.paramTrans(:,2) = 1./(1+exp(-params(:,2))); %learning rate       
        if m==4 %decay rate
            cbm.output.paramTrans(:,4) = 1./(1+exp(-params(:,4)));
        end
    
        for s=1:n_all
            subNb = subID_list(s);
            P = table2array(data(data.subNb==subNb,2:end));
            ntg = sum(P(:,19)==0);
            [ll,P_pred] = func_list_EL{m}(params(s,:),P);
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
    end
    
    fitRecapEL.pseudoR2{:,m} = cbm.output.pseudoR2;
    fitRecapEL.AIC{:,m}      = cbm.output.AIC;
    fitRecapEL.BIC{:,m}      = cbm.output.BIC;
    fitRecapEL.LHsub{:,m}    = cbm.output.LHsub;
    fitRecapEL.corrsub{:,m}  = cbm.output.corrsub;
    
    fitRecapEL.paramRaw.(mod_names_EL{m}) = cbm.output.parameters;
    fitRecapEL.paramTrans.(mod_names_EL{m}) = cbm.output.paramTrans; 
end
save([out_dir fs 'Recap_model_fitting_EL.mat'],'fitRecapEL');

%hierarchical model fitting
fname_hbi_EL = {'hbi_ExpLearn.mat';'hbi_ExpLearn_nomag.mat';'hbi_ExpLearn_mag.mat';'hbi_ExpLearn_decay.mat'};
parfor (m=2:n_mod_EL,3) 
    cbm_hbi(data_all, func_list_EL(m), {[out_dir fs out_fname_list_EL{m}]}, [out_dir fs fname_hbi_EL{m}]);
end

%calculate AIC and other measures based on hierarchical fits
hfitRecapEL = struct();
hfitRecapEL.pseudoR2 = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
hfitRecapEL.AIC = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
hfitRecapEL.BIC = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
hfitRecapEL.LHsub = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
hfitRecapEL.corrsub = table('Size',[n_all n_mod_EL],'VariableTypes',repmat({'double'},1,n_mod_EL),'VariableNames',mod_names_EL);
for m=1:n_mod_EL
    fname = [out_dir fs fname_hbi_EL{m}];
    load(fname,'cbm')
    params = cbm.output.parameters{:};
    npar = np(m);
    
    if m >= 2
        cbm.output.paramTrans = params;
        cbm.output.paramTrans(:,1) = exp(params(:,1)); %softmax beta
        cbm.output.paramTrans(:,2) = 1./(1+exp(-params(:,2))); %learning rate       
        if m==4 %decay rate
            cbm.output.paramTrans(:,4) = 1./(1+exp(-params(:,4)));
        end
    
        for s=1:n_all
            subNb = subID_list(s);
            P = table2array(data(data.subNb==subNb,2:end));
            ntg = sum(P(:,19)==0);
            [ll,P_pred] = func_list_EL{m}(params(s,:),P);
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
    end

    hfitRecapEL.pseudoR2{:,m} = cbm.output.pseudoR2;
    hfitRecapEL.AIC{:,m}      = cbm.output.AIC;
    hfitRecapEL.BIC{:,m}      = cbm.output.BIC;
    hfitRecapEL.LHsub{:,m}    = cbm.output.LHsub;
    hfitRecapEL.corrsub{:,m}  = cbm.output.corrsub;
    
    hfitRecapEL.paramRaw.(mod_names_EL{m}) = cbm.output.parameters;
    hfitRecapEL.paramTrans.(mod_names_EL{m}) = cbm.output.paramTrans;
end
save([out_dir fs 'Recap_model_fitting_EL.mat'],'fitRecapEL','hfitRecapEL');
%AIC values reported in Table S3:
AIC_values_EL = mean(table2array(fitRecapEL.AIC))