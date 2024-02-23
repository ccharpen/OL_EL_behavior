%This script performs the out-of-sample accuracy analyses for Study 2, and 
%generates the results presented in Tables 1 (OOS acc),
%as well as the plots presented in Figure 5B.

clear all
close all
fs = filesep;

data_dir = ['..' fs '..' fs 'data'];
mod_dir = ['..' fs 'dependencies' fs 'model_functions']; %directory where modelling functions are saved (common to both studies)

addpath(['..' fs 'dependencies']);
addpath(['..' fs 'dependencies' fs 'cbm-master' fs 'codes']);
addpath(['..' fs 'dependencies' fs 'plotSpread']);
addpath(mod_dir);

% init the randomization screen
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

%load data and format data as needed (in cell array)
data = readtable([data_dir fs 'data_study2.csv']);
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


%% Out-of-sample predictive accuracy
%calculate this between subjects across the entire sample AND within groups
%do this 6-fold

mod_names = {'Baseline'; 'ExpLearn'; 'ObsLearn'; 'FixArb'; 'DynArb'};
ll_mod_list = {@LL_Baseline; @LL_ExpLearn; @LL_ObsLearn; @LL_FixArb; @LL_DynArb};

%do one permutation for entire group first
id_list = unique(data.subNb);
acc_all = nan(100,5);
for p=1:100 %100 permutations
    p
    ord = [randperm(n_all)' id_list];
    ord = sortrows(ord,1);
    if rem(n_all,6) == 0
        ord(:,3) = repmat(1:6,1,floor(n_all/6))';
    else
        ord(:,3) = [repmat(1:6,1,floor(n_all/6))'; (1:rem(n_all,6))'];
    end
    ord = sortrows(ord,2);

    acc_perm = nan(6,5);
    for f=1:6
        train_id = ord(:,3)~=f;
        test_id = ord(ord(:,3)==f,2);

        acc = nan(length(test_id),5);
        for m=1:5
            params = fitRecap.paramRaw.(mod_names{m});
            params_train = mean(params(train_id,:));
            for s=1:length(test_id)
                sub = test_id(s);
                P = table2array(data(data.subNb==sub,2:end));
                [~, vals] = ll_mod_list{m}(params_train,P);
                acc(s,m) = nanmean(vals(:,1));
            end
        end
        acc_perm(f,:) = mean(acc);
    end
    acc_all(p,:) = mean(acc_perm);
end
OOSaccuracy.AllSubs = acc_all;

%group-by-group
for g=1:5
    g
    id_list_group = id_list(group == g);
    n_grp = length(id_list_group);
    
    acc_all = nan(100,5);
    parfor (p=1:100,5) %100 permutations
        ord = [randperm(n_grp)' id_list_group];
        ord = sortrows(ord,1);
        if rem(n_grp,6) == 0
            ord(:,3) = repmat(1:6,1,floor(n_grp/6))';
        else
            ord(:,3) = [repmat(1:6,1,floor(n_grp/6))'; (1:rem(n_grp,6))'];
        end
        ord = sortrows(ord,2);

        acc_perm = nan(6,5);
        for f=1:6
            train_id = ord(:,3)~=f;
            test_id = ord(ord(:,3)==f,2);

            acc = nan(length(test_id),5);
            for m=1:5
                params = fitRecap.paramRaw.(mod_names{m});
                params_train = mean(params(train_id,:));
                for s=1:length(test_id)
                    sub = test_id(s);
                    P = table2array(data(data.subNb==sub,2:end));
                    [~, vals] = ll_mod_list{m}(params_train,P);
                    acc(s,m) = nanmean(vals(:,1));
                end
            end
            acc_perm(f,:) = mean(acc);
        end
        acc_all(p,:) = mean(acc_perm);
    end
    
    if g==1
        OOSaccuracy.BaselineGroup = acc_all;
    elseif g==2
        OOSaccuracy.ExpLearnGroup = acc_all;
    elseif g==3
        OOSaccuracy.ObsLearnGroup = acc_all;
    elseif g==4
        OOSaccuracy.FixArbGroup = acc_all;
    elseif g==5
        OOSaccuracy.DynArbGroup = acc_all;
    end
    
end
save('Out_of_sample_accuracy.mat','OOSaccuracy')

recap = [mean(OOSaccuracy.AllSubs); mean(OOSaccuracy.BaselineGroup); ...
    mean(OOSaccuracy.ExpLearnGroup); mean(OOSaccuracy.ObsLearnGroup); ...
    mean(OOSaccuracy.FixArbGroup); mean(OOSaccuracy.DynArbGroup)];

figure;
h = heatmap(mod_names, [{'AllSubs'}; mod_names], recap);
h.Title = 'Between-subjects out of sample (6-fold) accuracy';
h.XLabel = 'Model';
h.YLabel = 'Group';
h.Colormap = redwhiteblue(0,1);


%% Within-subject predictive accuracy
%because of the clear breaks between blocks (i.e. new box stimuli), we can
%also do within subject OOS across blocks

%specify loglikelihood functions
func_list = {@LL_Baseline; @LL_ExpLearn; @LL_ObsLearn; @LL_FixArb; @LL_DynArb};
n_mod = length(func_list);

%specify output files
out_dir = 'oos_accuracy_outputs';
if ~exist(out_dir, 'dir')
   mkdir(out_dir)
end

mod_names = {'Baseline'; 'ExpLearn'; 'ObsLearn'; 'FixArb'; 'DynArb'};

%number of parameters for each model
np = [4;3;2;6;6]; 

%specify parameter priors
v = 6.25; %parameter variance (6.25 is large enough to cover a wide range of parameters with no excessive penalty)

for b=1:8
    %format data with subset of trials for training blocks
    data_all = cell(n_all,1);
    for s=1:n_all
        subNb = subID_list(s);
        data_all{s} = table2array(data(data.subNb==subNb & data.blockNb~=b,2:end));
    end
    
    %model fitting on training data
    parfor (m=1:n_mod,n_mod) 
        prior = struct('mean',zeros(np(m),1),'variance',v);
        fname = ['lap_' mod_names{m} '_testBlock' num2str(b) '.mat']
        cbm_lap(data_all, func_list{m}, prior, [out_dir fs fname]);
    end
end

acc = nan(n_all,n_mod,8);
for b=1:8
    %extract mean parameter and accuracy on test data
    for m=1:n_mod
        fname = [out_dir fs 'lap_' mod_names{m} '_testBlock' num2str(b) '.mat'];
        load(fname,'cbm')
        params = cbm.output.parameters;
        for s=1:n_all
            subNb = subID_list(s);
            P_test = table2array(data(data.subNb==subNb & data.blockNb==b,2:end));
            [~, vals] = ll_mod_list{m}(params(s,:),P_test);
            acc(s,m,b) = nanmean(vals(:,1));
        end
    end
end
OOSaccuracy.WithinSubs = mean(acc,3);
save('Out_of_sample_accuracy.mat','OOSaccuracy')

%% Plot Figure 5B
load('Out_of_sample_accuracy.mat','OOSaccuracy')

recap = [mean(OOSaccuracy.WithinSubs); mean(OOSaccuracy.WithinSubs(group==1,:)); ...
    mean(OOSaccuracy.WithinSubs(group==2,:)); mean(OOSaccuracy.WithinSubs(group==3,:)); ...
    mean(OOSaccuracy.WithinSubs(group==4,:)); mean(OOSaccuracy.WithinSubs(group==5,:))];

figure;
h = heatmap(mod_names, [{'AllSubs'}; mod_names], recap);
h.Title = 'Out-of-sample accuracy (8-fold, cross-blocks)';
h.XLabel = 'Model';
h.YLabel = 'Group';
h.Colormap = redwhiteblue(0,1);

%the values of OOS acc reported in Table 2 are contained in the first row
%of recap (mean(OOSaccuracy.WithinSubs))
