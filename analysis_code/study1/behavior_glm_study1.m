clear all
close all

fs = filesep;

dir_data = ['..' fs '..' fs 'data'];
addpath(['..' fs 'dependencies' fs 'plotSpread'])

% init the randomization screen
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

%load data
data = readtable([dir_data fs 'data_study1.csv']);

nsub = length(unique(data.subNb)); %number of subjects
ntr = 160; %number of trials

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
    %description of P:
    %column 1: block number
    %column 2: condition (1: TP 0.8, Mag LU; 2: TP 0.8, Mag HU; 3: TP 0.6, Mag LU; 4: TP 0.8, Mag HU)
    %column 3: trial number
    %column 4: trial number per block
    %column 5: common transition probability (associated with Box 1; 1:common, 0:rare)
    %column 6: whether common (1) or rare (0) transition happens
    %column 7: reward probability of orange token (1-blue token)
    %column 8: goal token (1: orange, 2: blue)
    %column 9: box chosen by partner (1: box 1 majority orange, 2: box 2 majority blue)
    %column 10: whether partner's choice is correct (1) or not (0)
    %column 11: obs token shown (1: orange, 2: blue)
    %column 12: whether token is goal (1) or not (0)
    %column 13: participant choice (1: left, 0: right)
    %column 14: participant choice (1: orange token, 0: blue token)
    %column 15: rt
    %column 16: whether participant choice is correct (1) or not (0)
    %column 17: token shown (1: orange, 2: blue)
    %column 18: outcome
    %column 19: missed trial (1)
    %column 20: action chosen by partner (1: left, 0: right)
    %column 21: reward if orange
    %column 22: reward magnitude

    %Number of good trials, Accuracy & RT - total & per condition
    i_good = P(:,19)==0;
    Behavior.Accuracy(sub,1) = mean(P(i_good,16));
    Behavior.RT(sub,1) = mean(P(i_good,15));
    Behavior.NMiss(sub,1) = sum(P(:,19)==1);
    
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
    
    %Number of good trials, Accuracy & RT - total & per trial-type
    Behavior.NumTr_cond(sub,:) = [sum(i_good & OL_lu & EL_lu & RM_h) sum(i_good & OL_lu & EL_lu & ~RM_h) ...
        sum(i_good & OL_lu & ~EL_lu & RM_h) sum(i_good & OL_lu & ~EL_lu & ~RM_h) ...
        sum(i_good & ~OL_lu & EL_lu & RM_h) sum(i_good & ~OL_lu & EL_lu & ~RM_h) ...
        sum(i_good & ~OL_lu & ~EL_lu & RM_h) sum(i_good & ~OL_lu & ~EL_lu & ~RM_h)];
    Behavior.Accuracy_cond(sub,:) = [mean(P(i_good & OL_lu & EL_lu & RM_h,16)) mean(P(i_good & OL_lu & EL_lu & ~RM_h,16)) ...
        mean(P(i_good & OL_lu & ~EL_lu & RM_h,16)) mean(P(i_good & OL_lu & ~EL_lu & ~RM_h,16)) ...
        mean(P(i_good & ~OL_lu & EL_lu & RM_h,16)) mean(P(i_good & ~OL_lu & EL_lu & ~RM_h,16)) ...
        mean(P(i_good & ~OL_lu & ~EL_lu & RM_h,16)) mean(P(i_good & ~OL_lu & ~EL_lu & ~RM_h,16))];
    Behavior.RT_cond(sub,:) = [mean(P(i_good & OL_lu & EL_lu & RM_h,15)) mean(P(i_good & OL_lu & EL_lu & ~RM_h,15)) ...
        mean(P(i_good & OL_lu & ~EL_lu & RM_h,15)) mean(P(i_good & OL_lu & ~EL_lu & ~RM_h,15)) ...
        mean(P(i_good & ~OL_lu & EL_lu & RM_h,15)) mean(P(i_good & ~OL_lu & EL_lu & ~RM_h,15)) ...
        mean(P(i_good & ~OL_lu & ~EL_lu & RM_h,15)) mean(P(i_good & ~OL_lu & ~EL_lu & ~RM_h,15))];
    
    % Check learning behavior by averaging accuracy for trial 1 to 8 after
    %a switch (or from beginning of block)
    tr_ind2 = zeros(ntr,1);
    tr_ind2(1:6) = (1:6)';
    for t=7:ntr
        if P(t,4)==1 || P(t,7)~=P(t-1,7)
            tr_ind2(t)=1;
        elseif P(t,4)==2 || tr_ind2(t-1)==1
            tr_ind2(t)=2;
        elseif P(t,4)==3 || tr_ind2(t-1)==2
            tr_ind2(t)=3;
        elseif P(t,4)==4 || tr_ind2(t-1)==3
            tr_ind2(t)=4;
        elseif P(t,4)==5 || tr_ind2(t-1)==4
            tr_ind2(t)=5;
        elseif P(t,4)==6 || tr_ind2(t-1)==5
            tr_ind2(t)=6;
        elseif P(t,4)==7 || tr_ind2(t-1)==6
            tr_ind2(t)=7;
        elseif P(t,4)==8 || tr_ind2(t-1)==7
            tr_ind2(t)=8;
        end
    end
    Behavior.Learning(sub,:) = [mean(P(i_good & tr_ind2==1,16)) mean(P(i_good & tr_ind2==2,16)) ...
        mean(P(i_good & tr_ind2==3,16)) mean(P(i_good & tr_ind2==4,16)) mean(P(i_good & tr_ind2==5,16)) ...
        mean(P(i_good & tr_ind2==6,16)) mean(P(i_good & tr_ind2==7,16)) mean(P(i_good & tr_ind2==8,16))];
    

    %% Calculate variables for the GLM and for behavioral signatures of the two strategies
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
        
    %define trials consistent with OL and with EL    
    ch_OL = (past_act==1 & past_tok==-1 & P(:,14)==1) | (past_act==-1 & past_tok==1 & P(:,14)==0);
    ch_EL = (past_act==1 & past_tok==-1 & P(:,14)==0) | (past_act==-1 & past_tok==1 & P(:,14)==1);
    
    %trials where OL and EL make different predictions
    ind_diff = ch_OL | ch_EL;
    Behavior.Prop_OL_EL_diff(sub,1) = mean(ind_diff);
    
    %% Propensity to choose according to OL (vs EL) overall and per condition
    Behavior.Prop_OL_ch(sub,:) = [mean(ch_OL(ind_diff)) sum(ind_diff)];
    
    %main effects
    Behavior.Prop_OL_ch_me(sub,:) = [mean(ch_OL(ind_diff & OL_lu)) mean(ch_OL(ind_diff & ~OL_lu)) ...
        mean(ch_OL(ind_diff & EL_lu)) mean(ch_OL(ind_diff & ~EL_lu)) ...
        mean(ch_OL(ind_diff & RM_h)) mean(ch_OL(ind_diff & ~RM_h)) ...
        sum(ind_diff & OL_lu) sum(ind_diff & ~OL_lu) sum(ind_diff & EL_lu) ...
        sum(ind_diff & ~EL_lu) sum(ind_diff & RM_h) sum(ind_diff & ~RM_h)];
    
    %2*2 breakdow, focusing on OL and EL uncertainty 
    Behavior.Prop_OL_ch_unc(sub,:) = [mean(ch_OL(ind_diff & OL_lu & EL_lu))  mean(ch_OL(ind_diff & OL_lu & ~EL_lu)) ...
        mean(ch_OL(ind_diff & ~OL_lu & EL_lu)) mean(ch_OL(ind_diff & ~OL_lu & ~EL_lu))];
    
    %2*2*2 breakdown (including magnitude)
    Behavior.Prop_OL_ch_cond(sub,:) = [mean(ch_OL(ind_diff & OL_lu & EL_lu & RM_h)) mean(ch_OL(ind_diff & OL_lu & EL_lu & ~RM_h)) ...
        mean(ch_OL(ind_diff & OL_lu & ~EL_lu & RM_h)) mean(ch_OL(ind_diff & OL_lu & ~EL_lu & ~RM_h)) ...
        mean(ch_OL(ind_diff & ~OL_lu & EL_lu & RM_h)) mean(ch_OL(ind_diff & ~OL_lu & EL_lu & ~RM_h)) ...
        mean(ch_OL(ind_diff & ~OL_lu & ~EL_lu & RM_h)) mean(ch_OL(ind_diff & ~OL_lu & ~EL_lu & ~RM_h))];
    Behavior.Ntr_OL_ch_cond(sub,:) = [sum(ind_diff & OL_lu & EL_lu & RM_h) sum(ind_diff & OL_lu & EL_lu & ~RM_h) ...
        sum(ind_diff & OL_lu & ~EL_lu & RM_h) sum(ind_diff & OL_lu & ~EL_lu & ~RM_h) ...
        sum(ind_diff & ~OL_lu & EL_lu & RM_h) sum(ind_diff & ~OL_lu & EL_lu & ~RM_h) ...
        sum(ind_diff & ~OL_lu & ~EL_lu & RM_h) sum(ind_diff & ~OL_lu & ~EL_lu & ~RM_h)];
    
    %effect of magnitude variance condition
    ilum = P(:,2)==1 | P(:,2)==3;
    Behavior.Prop_OL_ch_magvar(sub,:) = [mean(ch_OL(ind_diff & ilum)) mean(ch_OL(ind_diff & ~ilum))];
    
    %% build matrices for glms
    good = ~isnan(past_tok) & ~isnan(P(:,14)); 

    %glm matrix for assessing hybrid behavior
    glm_mat1 = [zscore(past_act(good)) zscore(past_tok_comb(good)) P(good,14)];
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
        [P(good & OL_lu,14);P(good & ~OL_lu,14)]];
    
    glm_mat_ELunc = [glm_mat_ELunc; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
        [zeros(sum(good & EL_lu),1); zscore(past_act(good & ~EL_lu))] ...
        [zscore(past_tok_comb(good & EL_lu)); zeros(sum(good & ~EL_lu),1)] ...
        [zeros(sum(good & EL_lu),1); zscore(past_tok_comb(good & ~EL_lu))] ...
        [P(good & EL_lu,14);P(good & ~EL_lu,14)]];
    
    glm_mat_Mag = [glm_mat_Mag; ones(sum(good),1)*subNb ...
        [zscore(past_act(good & RM_h)); zeros(sum(good & ~RM_h),1)] ...
        [zeros(sum(good & RM_h),1); zscore(past_act(good & ~RM_h))] ...
        [zscore(past_tok_comb(good & RM_h)); zeros(sum(good & ~RM_h),1)] ...
        [zeros(sum(good & RM_h),1); zscore(past_tok_comb(good & ~RM_h))] ...
        [P(good & RM_h,14);P(good & ~RM_h,14)]];

    %% measures of degenerate/baseline strategies
    Behavior.Porange(sub,1) = nanmean(P(:,14));
    Behavior.Pleft(sub,1) = nanmean(P(:,13));
    %calculate propensity to repeat own past action (left/right)
    repeatOwnAct = zeros(ntr-1,1);
    for t=1:ntr-1
        if P(t,13) == P(t+1,13)
            repeatOwnAct(t) = 1;
        elseif isnan(P(t+1,13))
            repeatOwnAct(t) = NaN;
        end
    end
    %calculate propensity to repeat partner's last action (left/right)
    repeatPartAct = zeros(ntr,1);
    for t=1:ntr        
        if P(t,20) == P(t,13)
            repeatPartAct(t) = 1;
        elseif isnan(P(t,13))
            repeatPartAct(t) = NaN;
        end
    end
    Behavior.PrepeatOwnAct(sub,1) = nanmean(repeatOwnAct);
    Behavior.PrepeatPartAct(sub,1) = nanmean(repeatPartAct);
    
    %% plot timeline of uncertainty conditions
    if sub==3
        figure;
        plot(P(:,7),'Color','#0072BD','LineWidth',2); hold on
        plot(~EL_lu*0.75+0.15,'.','Color','#0072BD','MarkerSize',8);
        plot(P(:,5),'Color','#D95319','LineWidth',2);
        plot(~OL_lu*0.85+0.1,'.','Color','#D95319','MarkerSize',8);
        plot(RM_h*0.95+0.05,'.','Color','#EDB120','MarkerSize',8);
        ylim([0 1])
        set(gca,'box','off')
        xlabel('Trial')
    end
    
end
save('Behavioral_variables.mat','Behavior')

%% Run mixed effect glms

%hybrid behavior
glm_tbl = array2table(glm_mat_allsubs, 'VariableNames',{'subNb','PastAct','PastTok','choice'});
glme = fitglme(glm_tbl,'choice ~ 1 + PastAct + PastTok + (1 + PastAct + PastTok|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
disp(glme)
anova(glme)
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
save('Behavioral_variables.mat','Behavior')

%separate GLMs for each uncertainty condition
%OL uncertainty
glm_OLunc_tbl = array2table(glm_mat_OLunc, 'VariableNames',...
    {'subNb','PastAct_LowOLU','PastAct_HighOLU','PastTok_LowOLU','PastTok_HighOLU','choice'});
glme_OLunc = fitglme(glm_OLunc_tbl,...
    'choice ~ 1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU + (1 + PastAct_LowOLU + PastAct_HighOLU + PastTok_LowOLU + PastTok_HighOLU|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
disp(glme_OLunc)
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
save('Behavioral_variables.mat','Behavior')

%EL uncertainty
glm_ELunc_tbl = array2table(glm_mat_ELunc, 'VariableNames',...
    {'subNb','PastAct_LowELU','PastAct_HighELU','PastTok_LowELU','PastTok_HighELU','choice'});
glme_ELunc = fitglme(glm_ELunc_tbl,...
    'choice ~ 1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU + (1 + PastAct_LowELU + PastAct_HighELU + PastTok_LowELU + PastTok_HighELU|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
disp(glme_ELunc)
[~,~,stats] = randomEffects(glme_ELunc);
Meffects = dataset2table(glme_ELunc.Coefficients);
Reffects = dataset2table(stats);
U = unstack(Reffects(:,2:4),'Estimate','Name');
U = renamevars(U,{'Level','x_Intercept_'},{'subNb','Intercept'});
U.Intercept = U.Intercept + Meffects.Estimate(1);
U.PastAct_LowELU = U.PastAct_LowELU + Meffects.Estimate(2);
U.PastAct_HighELU = U.PastAct_HighELU + Meffects.Estimate(3);
U.PastTok_LowELU = U.PastTok_LowELU + Meffects.Estimate(4);
U.PastTok_HighELU = U.PastTok_HighELU + Meffects.Estimate(5);
Behavior.GLME_ELunc.input_data = glm_ELunc_tbl;
Behavior.GLME_ELunc.Coefficients = Meffects;
Behavior.GLME_ELunc.Reffects = U;
Behavior.GLME_ELunc.R_stats = Reffects;
save('Behavioral_variables.mat','Behavior')

%reward magnitude
glm_Mag_tbl = array2table(glm_mat_Mag, 'VariableNames',...
    {'subNb','PastAct_HighMag','PastAct_LowMag','PastTok_HighMag','PastTok_LowMag','choice'});
glme_Mag = fitglme(glm_Mag_tbl,...
    'choice ~ 1 + PastAct_HighMag + PastAct_LowMag + PastTok_HighMag + PastTok_LowMag + (1 + PastAct_HighMag + PastAct_LowMag + PastTok_HighMag + PastTok_LowMag|subNb)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
disp(glme_Mag)
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
save('Behavioral_variables.mat','Behavior')

% association between glm effects and accuracy
ol_glme = Behavior.GLME.Reffects.PastAct;
el_glme = Behavior.GLME.Reffects.PastTok;
corr(ol_glme,Behavior.Accuracy)
corr(el_glme,Behavior.Accuracy)

%% Plots
%learning behavior
figure;  hold;
plot(1:8,mean(Behavior.Learning),'-k','LineWidth',2);
errorbar(1:8,mean(Behavior.Learning),std(Behavior.Learning)/sqrt(nsub),'.k','LineWidth',2);
xlim([0.5 8.5])
ylim([0.5 0.65])
xlabel('Trial since last switch')
ylabel('Mean Accuracy')
title('Learning behavior')

%main effects of condition on propensity to choose OL vs EL
data_OL_ch = [Behavior.Prop_OL_ch(:,1) Behavior.Prop_OL_ch_me(:,1:6)];
figure;
bar(1:7,mean(data_OL_ch),0.5,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1); hold on
plotSpread(data_OL_ch,'distributionColors',[0.4 0.4 0.4]);
errorbar(1:7,mean(data_OL_ch),std(data_OL_ch)/sqrt(nsub),'.k','LineWidth',1.5);
plot([0 8],[0.5 0.5],'--r')
xticklabels({'All trials','LU_o_b_s','HU_o_b_s','LU_e_x_p','HU_e_x_p','Hi_m_a_g','Lo_m_a_g'})
xlabel('Condition (trial-based)')
xtickangle(30)
ylabel('Proportion of OL choices (vs EL)')
ylim([0 1])

%2*2 breakdown
figure;
b = bar([nanmean(Behavior.Prop_OL_ch_unc(:,1:2)); nanmean(Behavior.Prop_OL_ch_unc(:,3:4))],0.8,...
    'FaceColor','flat','EdgeColor','k','LineWidth',1); hold on
b(1).CData = [0.9290 0.6940 0.1250]; 
b(2).CData = [0.4940 0.1840 0.5560];
plotSpread(Behavior.Prop_OL_ch_unc(:,1:4),'xValues',[0.85 1.15 1.85 2.15],'distributionColors',[0.4 0.4 0.4]); 
errorbar([0.85 1.15 1.85 2.15],nanmean(Behavior.Prop_OL_ch_unc(:,1:4)),...
    nanstd(Behavior.Prop_OL_ch_unc(:,1:4))/sqrt(nsub),'.k','LineWidth',1.5)
plot([0.5 2.5],[0.5 0.5],'--k')
xticks([1 2])
xticklabels({'Low','High'})
xlabel('OL uncertainty')
ylabel('Proportion of OL choices (vs EL)')
leg = legend({'Low','High'});
title(leg,'EL uncertainty')
set(gca,'box','off')

%2*2*2 breakdown
figure;
subplot(1,2,1); hold on
b = bar([nanmean(Behavior.Prop_OL_ch_cond(:,1:2)); nanmean(Behavior.Prop_OL_ch_cond(:,3:4))],0.8,...
    'FaceColor','flat','EdgeColor','k','LineWidth',1);
b(1).CData = [0.9290 0.6940 0.1250]; 
b(2).CData = [0.4940 0.1840 0.5560];
plotSpread(Behavior.Prop_OL_ch_cond(:,1:4),'xValues',[0.85 1.15 1.85 2.15],'distributionColors',[0.4 0.4 0.4]); 
errorbar([0.85 1.15 1.85 2.15],nanmean(Behavior.Prop_OL_ch_cond(:,1:4)),...
    nanstd(Behavior.Prop_OL_ch_cond(:,1:4))/sqrt(nsub),'.k','LineWidth',1.5)
plot([0.5 2.5],[0.5 0.5],'--k')
xticks([1 2])
xticklabels({'Low','High'})
xlabel('EL uncertainty')
ylabel('Proportion of OL choices (vs EL)')
title('LOW OL uncertainty')
subplot(1,2,2); hold on
b = bar([nanmean(Behavior.Prop_OL_ch_cond(:,5:6)); nanmean(Behavior.Prop_OL_ch_cond(:,7:8))],0.8,...
    'FaceColor','flat','EdgeColor','k','LineWidth',1);
b(1).CData = [0.9290 0.6940 0.1250]; 
b(2).CData = [0.4940 0.1840 0.5560];
plotSpread(Behavior.Prop_OL_ch_cond(:,5:8),'xValues',[0.85 1.15 1.85 2.15],'distributionColors',[0.4 0.4 0.4]); 
errorbar([0.85 1.15 1.85 2.15],nanmean(Behavior.Prop_OL_ch_cond(:,5:8)),...
    nanstd(Behavior.Prop_OL_ch_cond(:,5:8))/sqrt(nsub),'.k','LineWidth',1.5)
plot([0.5 2.5],[0.5 0.5],'--k')
xticks([1 2])
xticklabels({'Low','High'})
xlabel('EL uncertainty')
title('HIGH OL uncertainty')
leg = legend({'High','Low'});
title(leg,'Magnitude')

%GLM
figure; hold
bar(1,Behavior.GLME.Coefficients.Estimate([3 2]),0.57,'EdgeColor','k','LineWidth',1); 
plotSpread(Behavior.GLME.Reffects.PastTok,'xValues',0.85,'distributionColors',[0 0.32 0.47],'spreadWidth',0.25);
plotSpread(Behavior.GLME.Reffects.PastAct,'xValues',1.15,'distributionColors',[0.52 0.26 0.08],'spreadWidth',0.25);
errorbar([0.85 1.15],Behavior.GLME.Coefficients.Estimate([3 2]),Behavior.GLME.Coefficients.SE([3 2]),'.k','LineWidth',1.5);
xticks([])
ylabel('ME-GLM effect')
legend({'Effect of past outcome (EL)','Effect of past partner''s action (OL)'})

figure;
subplot(1,2,1); hold on
plot(Behavior.Prop_OL_ch(:,1),Behavior.GLME.Reffects.PastAct,'.'); lsline;
xlabel('Proportion of OL (vs EL) choices')
ylabel('Partner''s action ME-GLM effect')
title(['R = ' num2str(corr(Behavior.Prop_OL_ch(:,1),Behavior.GLME.Reffects.PastAct))])
subplot(1,2,2); hold on
plot(Behavior.Prop_OL_ch(:,1),Behavior.GLME.Reffects.PastTok,'.'); lsline;
xlabel('Proportion of OL (vs EL) choices')
ylabel('Previous outcome ME-GLM effect')
title(['R = ' num2str(corr(Behavior.Prop_OL_ch(:,1),Behavior.GLME.Reffects.PastTok))])

%GLM OL and EL uncertainty
figure;
subplot(1,2,1); hold on
h(1) = plot(1:2,Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),'-', 'LineWidth', 1.5);
h(2) = plot(1:2,Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),'-', 'LineWidth', 1.5);
plotSpread([Behavior.GLME_OLunc.Reffects.PastTok_LowOLU Behavior.GLME_OLunc.Reffects.PastTok_HighOLU],...
    'xValues',(1:2)-0.15,'distributionColors',[0 0.4470 0.7410],'spreadWidth',0.2);
plotSpread([Behavior.GLME_OLunc.Reffects.PastAct_LowOLU Behavior.GLME_OLunc.Reffects.PastAct_HighOLU],...
    'xValues',(1:2)+0.15,'distributionColors',[0.8500 0.3250 0.0980],'spreadWidth',0.2);
errorbar(1:2,Behavior.GLME_OLunc.Coefficients.Estimate([4 5]),Behavior.GLME_OLunc.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1:2,Behavior.GLME_OLunc.Coefficients.Estimate([2 3]),Behavior.GLME_OLunc.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35])
xticks([1 2])
xticklabels({'Low','High'})
xlabel('OL uncertainty')
ylabel('ME-GLM effect')
subplot(1,2,2); hold on
h(1) = plot(1:2,Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5);
h(2) = plot(1:2,Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_ELunc.Reffects.PastTok_LowELU Behavior.GLME_ELunc.Reffects.PastTok_HighELU],...
    'xValues',(1:2)-0.15,'distributionColors',[0 0.4470 0.7410],'spreadWidth',0.2);
plotSpread([Behavior.GLME_ELunc.Reffects.PastAct_LowELU Behavior.GLME_ELunc.Reffects.PastAct_HighELU],...
    'xValues',(1:2)+0.15,'distributionColors',[0.8500 0.3250 0.0980],'spreadWidth',0.2);
errorbar(1:2,Behavior.GLME_ELunc.Coefficients.Estimate([4 5]),Behavior.GLME_ELunc.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1:2,Behavior.GLME_ELunc.Coefficients.Estimate([2 3]),Behavior.GLME_ELunc.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35])
xticks([1 2])
xticklabels({'Low','High'})
xlabel('EL uncertainty')
ylabel('ME-GLM effect')
legend(h, {'Past outcome (EL effect)','Past action (OL effect)'});

%GLM reward magnitude 
figure;
h(1) = plot(1:2,Behavior.GLME_Mag.Coefficients.Estimate([4 5]),'-','LineWidth', 1.5); hold
h(2) = plot(1:2,Behavior.GLME_Mag.Coefficients.Estimate([2 3]),'-','LineWidth', 1.5);
plotSpread([Behavior.GLME_Mag.Reffects.PastTok_HighMag Behavior.GLME_Mag.Reffects.PastTok_LowMag],...
    'xValues',(1:2)-0.15,'distributionColors',[0 0.4470 0.7410],'spreadWidth',0.2);
plotSpread([Behavior.GLME_Mag.Reffects.PastAct_HighMag Behavior.GLME_Mag.Reffects.PastAct_LowMag],...
    'xValues',(1:2)+0.15,'distributionColors',[0.8500 0.3250 0.0980],'spreadWidth',0.2);
errorbar(1:2,Behavior.GLME_Mag.Coefficients.Estimate([4 5]),Behavior.GLME_Mag.Coefficients.SE([4 5]),...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5)
errorbar(1:2,Behavior.GLME_Mag.Coefficients.Estimate([2 3]),Behavior.GLME_Mag.Coefficients.SE([2 3]),...
    'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
xlim([0.60 2.35])
xticks([1 2])
xticklabels({'High','Low'})
xlabel('Reward Magnitude')
ylabel('ME-GLM effect')
legend(h, {'Past outcome (EL effect)','Past action (OL effect)'});
