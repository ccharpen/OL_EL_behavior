clear all
close all

fs = filesep;

addpath(['..' fs 'dependencies' fs 'plotSpread'])

%load behavioral variable and model-fitting outputs
load('Behavioral_variables.mat')
out_dir = 'model_fitting_outputs';
load([out_dir fs 'Recap_model_fitting.mat'])
load([out_dir fs 'hbi_5mods.mat'])

n_all = length(Behavior.subID_list);

%define groups based on hierarchical fit model responsibility values
hfit_recap = cbm.output.responsibility;
for s=1:n_all
    hfit_recap(s,6) = find(hfit_recap(s,1:5)==max(hfit_recap(s,1:5)));
end
group = hfit_recap(:,6);
gsize = [sum(group==1) sum(group==2) sum(group==3) sum(group==4) sum(group==5)];
gdist = gsize/sum(gsize);
save('Recap_model_fitting.mat','fitRecap','hfitRecap','group');

%color palette for plots
clr = [248/255 125/255 115/255; 184/255 186/255 65/255; ...
    51/255 198/255 142/255; 34/255 181/255 246/255; ...
    239/255 110/255 253/255]; 

%calculate variable of interests to examine across groups
arb_index = Behavior.Prop_OL_ch_cond(:,4) - Behavior.Prop_OL_ch_cond(:,5);
base_index = mean([abs(Behavior.Porange-0.5) abs(Behavior.Pleft-0.5) ...
    abs(Behavior.PrepeatOwnAct-0.5) abs(Behavior.PrepeatPartAct-0.5)],2);
mean_learning = nan(8,5); sem_learning = nan(8,5);
mean_glm = nan(2,5); sem_glm = nan(2,5);
mean_OLch = nan(1,5); sem_OLch = nan(1,5);
mean_OLch_cond = nan(8,5); sem_OLch_cond = nan(8,5);
mean_arb = nan(1,5); sem_arb = nan(1,5);
mean_base = nan(1,5); sem_base = nan(1,5);
for g=1:5
    mean_learning(:,g) = nanmean(Behavior.Learning(group==g,:))';
    sem_learning(:,g) = nanstd(Behavior.Learning(group==g,:))'/sqrt(gsize(g));
    
    mean_glm(:,g) = [nanmean(Behavior.GLME.Reffects.PastTok(group==g)); ...
        nanmean(Behavior.GLME.Reffects.PastAct(group==g))];
    sem_glm(:,g) = [nanstd(Behavior.GLME.Reffects.PastTok(group==g))/sqrt(gsize(g)); ...
        nanstd(Behavior.GLME.Reffects.PastAct(group==g))/sqrt(gsize(g))];
        
    mean_OLch(g) = nanmean(Behavior.Prop_OL_ch(group==g,1));
    sem_OLch(g) = nanstd(Behavior.Prop_OL_ch(group==g,1))/sqrt(gsize(g));
    
    mean_OLch_cond(:,g) = nanmean(Behavior.Prop_OL_ch_cond(group==g,:))';
    sem_OLch_cond(:,g) = nanstd(Behavior.Prop_OL_ch_cond(group==g,:))'/sqrt(gsize(g));
        
    mean_arb(g) = nanmean(arb_index(group==g));
    sem_arb(g) = nanstd(arb_index(group==g))/sqrt(gsize(g));

    mean_base(g) = nanmean(base_index(group==g));
    sem_base(g) = nanstd(base_index(group==g))/sqrt(gsize(g));
end

%plot learning curves per group
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
xlabel('Trial since last switch')
ylabel('Accuracy')
ylim([0.4 0.85])
legend({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})

%plot GLM main effects per group
figure;
subplot(1,2,1); hold on
bar(1:5,mean_glm(1,:),'EdgeColor','k','LineWidth',1);
plotSpread(Behavior.GLME.Reffects.PastTok,'distributionIdx',group,'xValues',(1:5),'distributionColors',[0 0.32 0.47]);
errorbar((1:5),mean_glm(1,:),sem_glm(1,:),'.k','LineWidth',1.5);
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
ylim([-0.5 2])
xlabel('Group')
ylabel('ME-GLM effect')
title('Effect of past outcome (EL)')
subplot(1,2,2); hold on
bar(1:5,mean_glm(2,:),'FaceColor','#D95319','EdgeColor','k','LineWidth',1);
plotSpread(Behavior.GLME.Reffects.PastAct,'distributionIdx',group,'xValues',(1:5),'distributionColors',[0.52 0.26 0.08]);
errorbar((1:5),mean_glm(2,:),sem_glm(2,:),'.k','LineWidth',1.5);
ylim([-0.3 1])
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('ME-GLM effect')
title('Effect of past partner''s action (OL)')

%propensity to choose according to OL vs EL
figure; hold
b = bar(1:5,mean_OLch,0.7,'facecolor', 'flat','EdgeColor','k','LineWidth',1);
b.CData = clr;
plotSpread(Behavior.Prop_OL_ch(:,1),'distributionIdx',group,'distributionColors','k');
errorbar(1:5,mean_OLch,sem_OLch,'.k','LineWidth',1.5);
plot([0 6], [0.5 0.5], 'k--')
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('OL (vs EL) choice propensity')
ylim([0 1])

%arbitration index per group
figure; hold
b = bar(1:5,mean_arb,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
b.CData = clr;
plotSpread(arb_index,'distributionIdx',group,'distributionColors','k');
errorbar(1:5,mean_arb,sem_arb,'.k','LineWidth',1.5);
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('Arbitration index')
ylim([-0.8 1.2])

%plot GLM effects by OL and EL uncertainty per group
recap_pastAct_OLU = [Behavior.GLME_OLunc.Reffects.PastAct_LowOLU Behavior.GLME_OLunc.Reffects.PastAct_HighOLU];
recap_pastAct_ELU = [Behavior.GLME_ELunc.Reffects.PastAct_LowELU Behavior.GLME_ELunc.Reffects.PastAct_HighELU];
recap_pastTok_OLU = [Behavior.GLME_OLunc.Reffects.PastTok_LowOLU Behavior.GLME_OLunc.Reffects.PastTok_HighOLU];
recap_pastTok_ELU = [Behavior.GLME_ELunc.Reffects.PastTok_LowELU Behavior.GLME_ELunc.Reffects.PastTok_HighELU];
data_cell = {recap_pastAct_OLU, recap_pastAct_ELU, recap_pastTok_OLU, recap_pastTok_ELU};
xlabel_list = {'OL uncertainty','EL uncertainty','OL uncertainty','EL uncertainty'};
figure;
for sp=1:4
    data_plot = data_cell{sp};
    subplot(2,2,sp); hold on
    for g=1:5
        h(g) = plot(1:2,mean(data_plot(group==g,:)),'-','Color',clr(g,:),'LineWidth',1);
        plot(2.5,mean(data_plot(group==g,1)-data_plot(group==g,2)),'.','Color',clr(g,:))
        errorbar(1:2,mean(data_plot(group==g,:)),std(data_plot(group==g,:))/sqrt(sum(group==g)),...
            'Color',clr(g,:),'LineWidth',1)
        errorbar(2.5,mean(data_plot(group==g,1)-data_plot(group==g,2)),...
            std(data_plot(group==g,1)-data_plot(group==g,2))/sqrt(sum(group==g)),...
            'Color',clr(g,:),'LineWidth',1)
    end
    plot([0.75 2.75],[0 0],'--k')
    xlim([0.75 2.75])
    xticks([1 2 2.5])
    xticklabels({'Low','High','Diff'})
    xlabel(xlabel_list{sp})
    if sp<=2
        if sp==1
            ylabel('ME-GLM effect of past action (OL)')
        end
        ylim([-0.5 1.5])
    elseif sp>=3
        if sp==3
            ylabel('ME-GLM effect of past outcome (EL)')
        end
        ylim([-0.2 2])
    end
    if sp==4
        leg = legend(h, {'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'});
        title(leg,'group')
    end
end

%plot GLM effects by reward magnitude and group
recap_pastAct_Mag = [Behavior.GLME_Mag.Reffects.PastAct_HighMag Behavior.GLME_Mag.Reffects.PastAct_LowMag];
recap_pastTok_Mag = [Behavior.GLME_Mag.Reffects.PastTok_HighMag Behavior.GLME_Mag.Reffects.PastTok_LowMag];
figure;
subplot(2,1,1); hold on
for g=1:5
    h(g) = plot(1:2,mean(recap_pastAct_Mag(group==g,:)),'-','Color',clr(g,:),'LineWidth',1);
    plot(2.5,mean(recap_pastAct_Mag(group==g,1)-recap_pastAct_Mag(group==g,2)),'.','Color',clr(g,:))
    errorbar(1:2,mean(recap_pastAct_Mag(group==g,:)),std(recap_pastAct_Mag(group==g,:))/sqrt(sum(group==g)),...
        'Color',clr(g,:),'LineWidth',1)
    errorbar(2.5,mean(recap_pastAct_Mag(group==g,1)-recap_pastAct_Mag(group==g,2)),...
        std(recap_pastAct_Mag(group==g,1)-recap_pastAct_Mag(group==g,2))/sqrt(sum(group==g)),...
        'Color',clr(g,:),'LineWidth',1)
end
plot([0.75 2.75],[0 0],'--k')
xlim([0.75 2.75])
xticks([1 2 2.5])
xticklabels({'High','Low','Diff'})
xlabel('Reward Magnitude')
ylabel('ME-GLM effect of past action (OL)')
ylim([-0.5 1.5])
subplot(2,1,2); hold on
for g=1:5
    h(g) = plot(1:2,mean(recap_pastTok_Mag(group==g,:)),'-','Color',clr(g,:),'LineWidth',1);
    plot(2.5,mean(recap_pastTok_Mag(group==g,1)-recap_pastTok_Mag(group==g,2)),'.','Color',clr(g,:))
    errorbar(1:2,mean(recap_pastTok_Mag(group==g,:)),std(recap_pastTok_Mag(group==g,:))/sqrt(sum(group==g)),...
        'Color',clr(g,:),'LineWidth',1)
    errorbar(2.5,mean(recap_pastTok_Mag(group==g,1)-recap_pastTok_Mag(group==g,2)),...
        std(recap_pastTok_Mag(group==g,1)-recap_pastTok_Mag(group==g,2))/sqrt(sum(group==g)),...
        'Color',clr(g,:),'LineWidth',1)
end
plot([0.75 2.75],[0 0],'--k')
xlim([0.75 2.75])
xticks([1 2 2.5])
xticklabels({'High','Low','Diff'})
xlabel('Reward Magnitude')
ylabel('ME-GLM effect of past outcome (EL)')
ylim([-0.2 2])
leg = legend(h, {'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'});
title(leg,'group')

%baseline index per group
figure; hold
b = bar(1:5,mean_base,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
b.CData = clr;
plotSpread(base_index,'distributionIdx',group,'distributionColors','k');
errorbar(1:5,mean_base,sem_base,'.k','LineWidth',1.5);
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('Baseline strategies index')

    
%% plot group difference in IQ, gender, age, education
%load recap table
recap_table = readtable('recap_analyses_study2.csv');

%mean age per group
mean_age = [mean(recap_table.age(recap_table.group==1)) mean(recap_table.age(recap_table.group==2)) ...
    mean(recap_table.age(recap_table.group==3)) mean(recap_table.age(recap_table.group==4)) ...
    mean(recap_table.age(recap_table.group==5))];
sem_age = [std(recap_table.age(recap_table.group==1))/sqrt(gsize(1)) std(recap_table.age(recap_table.group==2))/sqrt(gsize(2)) ...
    std(recap_table.age(recap_table.group==3))/sqrt(gsize(3)) std(recap_table.age(recap_table.group==4))/sqrt(gsize(3)) ...
    std(recap_table.age(recap_table.group==5))/sqrt(gsize(5))];
figure; hold
b = bar(1:5,mean_age,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
b.CData = clr;
plotSpread(recap_table.age,'distributionIdx',recap_table.group,'distributionColors','k');
errorbar(1:5,mean_age,sem_age,'.k','LineWidth',1.5);
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('Age')
ylim([15 70])

%gender distribution per group
recap_gender = nan(5,3);
for g=1:5
    recap_gender(g,1) = sum(strcmp(recap_table.gender(recap_table.group==g),'F'))/gsize(g);
    recap_gender(g,2) = sum(strcmp(recap_table.gender(recap_table.group==g),'N'))/gsize(g);
    recap_gender(g,3) = sum(strcmp(recap_table.gender(recap_table.group==g),'M'))/gsize(g);
end
figure;
b = bar(recap_gender,0.7,'stacked','FaceColor','flat');
b(1).CData = [0.3 0.3 0.3];
b(2).CData = [0.5 0.5 0.5];
b(3).CData = [0.7 0.7 0.7];
leg = legend({'Female','Non-binary','Male'});
title(leg,'gender');
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('Proportion')

%mean IQ score per group
mean_ICAR = [nanmean(recap_table.ICAR_score(recap_table.group==1)) nanmean(recap_table.ICAR_score(recap_table.group==2)) ...
    nanmean(recap_table.ICAR_score(recap_table.group==3)) nanmean(recap_table.ICAR_score(recap_table.group==4)) ...
    nanmean(recap_table.ICAR_score(recap_table.group==5))];
sem_ICAR = [nanstd(recap_table.ICAR_score(recap_table.group==1))/sqrt(gsize(1)) nanstd(recap_table.ICAR_score(recap_table.group==2))/sqrt(gsize(2)) ...
    nanstd(recap_table.ICAR_score(recap_table.group==3))/sqrt(gsize(3)) nanstd(recap_table.ICAR_score(recap_table.group==4))/sqrt(gsize(3)) ...
    nanstd(recap_table.ICAR_score(recap_table.group==5))/sqrt(gsize(5))];
figure; hold
b = bar(1:5,mean_ICAR,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
b.CData = clr;
plotSpread(recap_table.ICAR_score,'distributionIdx',recap_table.group,'distributionColors','k');
errorbar(1:5,mean_ICAR,sem_ICAR,'.k','LineWidth',1.5);
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('ICAR score (~IQ)')

%education level per group
%High school diploma/GED = 1
%Associate's degree (2 yr) = 2
%Currently in college/university = 3
%Bachelor's degree (4 yr) = 4
%Master's degree = 5
%Doctoral degree = 6
recap_education = nan(5,6);
for g=1:5
    for e=1:6
        recap_education(g,e) = sum(recap_table.education(recap_table.group==g)==e)/gsize(g);
    end
end
figure;
bar(recap_education,0.7,'stacked');
leg = legend({'High-School','Associate','In college','Bachelors','Masters','Doctorate'});
title(leg,'education');
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('Proportion')
ylim([0 1])

mean_ed = [nanmean(recap_table.education(recap_table.group==1)) nanmean(recap_table.education(recap_table.group==2)) ...
    nanmean(recap_table.education(recap_table.group==3)) nanmean(recap_table.education(recap_table.group==4)) ...
    nanmean(recap_table.education(recap_table.group==5))];
sem_ed = [nanstd(recap_table.education(recap_table.group==1))/sqrt(gsize(1)) nanstd(recap_table.education(recap_table.group==2))/sqrt(gsize(2)) ...
    nanstd(recap_table.education(recap_table.group==3))/sqrt(gsize(3)) nanstd(recap_table.education(recap_table.group==4))/sqrt(gsize(3)) ...
    nanstd(recap_table.education(recap_table.group==5))/sqrt(gsize(5))];
figure; hold
b = bar(1:5,mean_ed,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
b.CData = clr;
plotSpread(recap_table.education,'distributionIdx',recap_table.group,'distributionColors','k');
errorbar(1:5,mean_ed,sem_ed,'.k','LineWidth',1.5);
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
xlabel('Group')
ylabel('Education level')

%plot parameter estimates per group
for p=53:61
    if p<=56
        par = abs(table2array(recap_table(:,p)));
        par_name = ['abs(' recap_table.Properties.VariableNames{p} ')'];
    else
        par = table2array(recap_table(:,p));
        par_name = recap_table.Properties.VariableNames{p};
    end
    mean_par = [nanmean(par(recap_table.group==1)) nanmean(par(recap_table.group==2)) ...
        nanmean(par(recap_table.group==3)) nanmean(par(recap_table.group==4)) ...
        nanmean(par(recap_table.group==5))];
    sem_par = [nanstd(par(recap_table.group==1))/sqrt(gsize(1)) nanstd(par(recap_table.group==2))/sqrt(gsize(2)) ...
        nanstd(par(recap_table.group==3))/sqrt(gsize(3)) nanstd(par(recap_table.group==4))/sqrt(gsize(3)) ...
        nanstd(par(recap_table.group==5))/sqrt(gsize(5))];
    
    figure; hold
    b = bar(1:5,mean_par,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
    b.CData = clr;
    plotSpread(par,'distributionIdx',recap_table.group,'distributionColors','k');
    errorbar(1:5,mean_par,sem_par,'.k','LineWidth',1.5);
    xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
    xtickangle(30)
    xlabel('Group')
    ylabel(par_name)
end