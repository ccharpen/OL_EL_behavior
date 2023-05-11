clear all
close all

fs = filesep;

addpath(['dependencies' fs 'cbm-master' fs 'codes']);
addpath(['dependencies' fs 'model_functions'])

%load data from study 1 to get an example dataset from each trial list
data_dir = ['..' fs 'data'];
data = readtable([data_dir fs 'data_study1.csv']);
%select 1 trial list of each type from P matrix (from study 1 data summary, 
%16 different trial lists in total)
ind_trial_list = [10;27;8;6;18;16;5;2;12;4;25;29;19;9;20;21];

%select models
ll_mod_list = {@LL_Baseline; @LL_ExpLearn; @LL_ObsLearn; @LL_FixArb; @LL_DynArb};

gc_mod_list = {@generate_choice_Baseline; @generate_choice_ExpLearn; ...
    @generate_choice_ObsLearn; @generate_choice_FixArb; @generate_choice_DynArb};

out_names = {'Baseline'; 'ExpLearn'; 'ObsLearn'; 'FixArb'; 'DynArb'};

v = 6.25; %parameter variance (6.25 is large enough to cover a wide range of parameters with no excessive penalty)

n_mod = length(ll_mod_list);
np = [4;3;2;6;6]; 
nsim = 100;


%% extract existing trial list and generate data
tr_list = zeros(nsim,1);
for n = 1:nsim
    tr_list(n) = ind_trial_list(randi(16));
end

%sample parameters (use high beta values for all models)
beta_OL_list = randn(nsim,1)+2;
beta_EL_list = randn(nsim,1)+2;
params_list = {randn(nsim,np(1)); ...
    [beta_EL_list randn(nsim,np(2)-1)]; ...
    [beta_OL_list randn(nsim,np(3)-1)]; ...
    [beta_OL_list beta_EL_list randn(nsim,np(4)-2)]; ...
    [beta_OL_list beta_EL_list randn(nsim,np(5)-2)]};

Recap_recov_params = cell(n_mod, n_mod);
Recap_pseudoR2 = cell(n_mod, n_mod);
Recap_recov_params_hbi = cell(n_mod, 1);
Recap_pseudoR2_hbi = cell(n_mod, 1);

%generate data for each model/each simulated set of parameters
data_sim_all = cell(nsim,n_mod);
for m = 1:n_mod
    params_mod = params_list{m};
    for n = 1:nsim
        P = table2array(data(data.subNb == tr_list(n),2:end));
        param_sim = params_mod(n,:);
        P_pred = gc_mod_list{m}(param_sim,P);
        P_new = P;
        P_new(:,14) = P_pred(:,4); %choice (1: orange, 0: blue)
        P_new(:,17) = P_pred(:,5); %token  (1: orange, 2: blue)
        P_new(:,16) = P_pred(:,6); %whether predicted choice is correct (1) or not (0)
        P_new(:,18) = P_pred(:,8)*100; %predicted outcome, should be scaled from 0 to 100
        if m==1
            P_new(:,13) = P_pred(:,9); %action (1: left, 0: right)
            %not needed for other models than baseline
        end
        data_sim_all{n,m} = P_new;
    end
end

%% fit all simulated datasets with all models

%optimize number of cores to use for parallelizing
numcores = feature('numcores');
if numcores>=5
    final_nc = 5;
else
    final_nc = numcores;
end

%create result directory if it doesn't exist
if ~exist('sim_results', 'dir')
   mkdir('sim_results')
end

for ms = 1:n_mod
    data_sim_all_mod = data_sim_all(:,ms);
    parfor (mf = 1:n_mod, final_nc)
        npar = np(mf);
        prior = struct('mean',zeros(npar,1),'variance',v);
        out_fname = ['sim_results' fs 'sim' out_names{ms} '_mod' out_names{mf} '.mat'];
        cbm_lap(data_sim_all_mod, ll_mod_list{mf}, prior, out_fname);
        res = load(out_fname);
        recov_params = res.cbm.output.parameters;
        pseudoR2 = zeros(nsim,1);
        for n=1:nsim
            P = data_sim_all_mod{n};
            [ll,~] = ll_mod_list{mf}(recov_params(n,:),P);
            pseudoR2(n) = 1 - ll/(160*log(0.5));
            if pseudoR2(n)<0
                pseudoR2(n)=0;
            end
        end
        Recap_recov_params{ms,mf} = recov_params;
        Recap_pseudoR2{ms,mf} = pseudoR2;
    end
end
%when simulated and recovered model are the same, run hbi on single model for better parameter estimates
nit=10;
parfor (m=1:n_mod, 4)
    npar = np(m);
    data_sim_all_mod = data_sim_all(:,m);
    out_fname = ['sim_results' fs 'sim' out_names{m} '_mod' out_names{m} '.mat'];
    recov_params = zeros(nsim,npar,nit);
    pseudoR2 = zeros(nsim,nit);
    for it=1:nit
        hbi_fname = ['sim_results' fs 'hbi_sim' out_names{m} '_mod' out_names{m} '_it' num2str(it) '.mat'];
        cbm_hbi(data_sim_all_mod, ll_mod_list(m), {out_fname}, {hbi_fname});
        res = load(hbi_fname);
        recov_params(:,:,it) = res.cbm.output.parameters{1,1};
        for n=1:nsim
            P = data_sim_all_mod{n};
            [ll,~] = ll_mod_list{m}(recov_params(n,:,it),P);
            pseudoR2(n,it) = 1 - ll/(160*log(0.5));
            if pseudoR2(n,it)<0
                pseudoR2(n,it)=0;
            end
        end
    end
    Recap_recov_params_hbi{m} = mean(recov_params,3);
    Recap_pseudoR2_hbi{m} = mean(pseudoR2,2);
end
Actual_params = params_list;
save(['sim_results' fs 'Recap_params_R2.mat'],'Actual_params','Recap_recov_params','Recap_pseudoR2',...
    'Recap_recov_params_hbi','Recap_pseudoR2_hbi','tr_list','data_sim_all')

%% run hierarchical fitting to plot exceedance probability
models = {@LL_Baseline; @LL_ExpLearn; @LL_ObsLearn; @LL_FixArb; @LL_DynArb};
mod_names = {'Baseline'; 'ExpLearn'; 'ObsLearn'; 'FixArb'; 'DynArb'};
cd sim_results
parfor (r=1:5, 3)
    data_sim_all_mod = data_sim_all(:,r);
    fcbm_maps = {['sim' mod_names{r} '_modBaseline.mat'];['sim' mod_names{r} '_modExpLearn.mat'];...
        ['sim' mod_names{r} '_modObsLearn.mat'];['sim' mod_names{r} '_modFixArb.mat'];...
        ['sim' mod_names{r} '_modDynArb.mat']};
    fname_hbi = ['hbi_sim' mod_names{r} '.mat'];
    cbm_hbi(data_sim_all_mod, models, fcbm_maps, fname_hbi);
end

EP_matrix = zeros(n_mod,n_mod);
for ms = 1:n_mod
    load(['hbi_sim' mod_names{ms} '.mat'])
    EP_matrix(ms,:) = cbm.output.exceedance_prob;
end
figure;
h = heatmap(mod_names,mod_names,EP_matrix);
h.Title = 'Exceedance Probability';
h.XLabel = 'Model used for fitting';
h.YLabel = 'Model used for data generation';
