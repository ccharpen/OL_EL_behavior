clear all
close all

fs = filesep;
addpath(['..' fs 'dependencies' fs 'plotSpread']);

%load data from each study and pool
recap_s1 = readtable(['..' fs 'study1' fs 'recap_analyses_study1.csv']);
recap_s2 = readtable(['..' fs 'study2' fs 'recap_analyses_study2.csv']);
recap_all = [recap_s1;recap_s2];
recap_all.study = [ones(height(recap_s1),1);2*ones(height(recap_s2),1)];
n_all = height(recap_all);
recap_all.uniqueID = (1:n_all)';

%load factor scores
factors = readtable('factor_scores_excl_8F.csv');

%load items
items = readtable(['..' fs '..' fs 'data' fs 'total_questionnaire_items.csv']);

recap_all.F1_depression = nan(n_all,1);
recap_all.F2_socialAnx = nan(n_all,1);
recap_all.F3_autism = nan(n_all,1);
recap_all.F4_stateAnx = nan(n_all,1);
recap_all.F5_socialResp = nan(n_all,1);
recap_all.F6_groupAvoid = nan(n_all,1);
recap_all.F7_traitAnx = nan(n_all,1);
recap_all.F8_perfAnx = nan(n_all,1);
for s=1:n_all
    sID = recap_all.uniqueID(s);
    inds = find(factors.uniqueID==sID);
    inds_item = find(items.uniqueID==sID);
    if ~isempty(inds)
        recap_all.F1_depression(s) = factors.WLS1(inds);
        recap_all.F2_socialAnx(s) = factors.WLS2(inds);
        recap_all.F3_autism(s) = factors.WLS3(inds);
        recap_all.F4_stateAnx(s) = factors.WLS4(inds);
        recap_all.F5_socialResp(s) = factors.WLS5(inds);
        recap_all.F6_groupAvoid(s) = factors.WLS6(inds);
        recap_all.F7_traitAnx(s) = factors.WLS7(inds);
        recap_all.F8_perfAnx(s) = factors.WLS8(inds);
    end
end
excl = isnan(recap_all.F1_depression);
recap_excl = recap_all(~excl,:);
recap_excl.zICAR = zscore(recap_excl.ICAR_score); %add zscore of ICAR
items_excl = items(~excl,:);
recap_excl = [recap_excl(:,80) recap_excl(:,1) recap_excl(:,81) recap_excl(:,2:6) recap_excl(:,90) recap_excl(:,7:79) recap_excl(:,82:89)]; %reorder variables
writetable(recap_excl,'recap_factor_groups.csv')

%colors for bar plots
clr = [248/255 125/255 115/255; 184/255 186/255 65/255; ...
    51/255 198/255 142/255; 34/255 181/255 246/255; ...
    239/255 110/255 253/255];

%group sizes
grp = recap_excl.group;
gsize = [sum(grp==1) sum(grp==2) sum(grp==3) sum(grp==4) sum(grp==5)];

%plot factor profile for each group
g_list = {'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'};
figure;
for g = 1:5
    vals = table2array(recap_excl(grp==g,83:90));
    subplot(2,3,g); hold on
    b = bar(1:8,mean(vals),0.7,'facecolor',clr(g,:),'EdgeColor','k','LineWidth',1);
    %plotSpread(vals,'distributionIdx',grp,'distributionColors','k');
    errorbar(1:8,mean(vals),std(vals)/sqrt(gsize(g)),'.k','LineWidth',1.5);
    xticks(1:8)
    xticklabels({'depression','socialAnx','autism','stateAnx','socialResp','groupAvoid','traitAnx','perfAnx'})
    xtickangle(30)
    ylabel('factor score')
    title([g_list{g} ' group'])
end
subplot(2,3,6); hold on
for g=1:5
    x = mean(recap_excl.F3_autism(grp==g));
    y = mean(recap_excl.F7_traitAnx(grp==g));
    xerr = std(recap_excl.F3_autism(grp==g))/sqrt(gsize(g));
    yerr = std(recap_excl.F7_traitAnx(grp==g))/sqrt(gsize(g));
    plot(x,y,'.','Color',clr(g,:),'MarkerSize',20)
end
for g=1:5
    x = mean(recap_excl.F3_autism(grp==g));
    y = mean(recap_excl.F7_traitAnx(grp==g));
    xerr = std(recap_excl.F3_autism(grp==g))/sqrt(gsize(g));
    yerr = std(recap_excl.F7_traitAnx(grp==g))/sqrt(gsize(g));
    errorbar(x,y,yerr,yerr,xerr,xerr,'Color',clr(g,:))
end
xlabel("mean autism factor score");
ylabel("mean traitAnx factor score");
h = legend({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'});
title(h,'group')

%alternative visualization: plot group means for each factor
figure;
for sp = 1:8
    tcol = 82+sp;
    vname = recap_excl.Properties.VariableNames{tcol};
    vals = recap_excl.(vname);
    mean_f = [mean(vals(grp==1)) mean(vals(grp==2)) mean(vals(grp==3)) mean(vals(grp==4)) mean(vals(grp==5))];
    sem_f = [std(vals(grp==1))/sqrt(gsize(1)) std(vals(grp==2))/sqrt(gsize(2)) ...
        std(vals(grp==3))/sqrt(gsize(3)) std(vals(grp==4))/sqrt(gsize(4)) std(vals(grp==5))/sqrt(gsize(5))];
    subplot(2,4,sp); hold on
    b = bar(1:5,mean_f,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
    b.CData = clr;
    %plotSpread(vals,'distributionIdx',grp,'distributionColors','k');
    errorbar(1:5,mean_f,sem_f,'.k','LineWidth',1.5);
    xticks(1:5)
    xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
    xtickangle(30)
    title([vname(1:2) ' - ' vname(4:end)])
    if sp>=5
        xlabel('Group')
    end
end


%% plot F3_autism and F7_traitAnxiety per group
figure;
x=recap_excl.F3_autism;
y=recap_excl.F7_traitAnx;
[p,S] = polyfit(x,y,1);
xv = linspace(min(x), max(x), 150);
[y_ext,delta] = polyconf(p,xv,S,'predopt','curve');
scatterhist(x,y,'Group',recap_excl.group,'Kernel','on','LineStyle','-','Color',clr,'Marker','.')
hold on
plot(xv, y_ext, '-k', 'LineWidth', 1.2)
patch([xv fliplr(xv)], [(y_ext+delta) fliplr((y_ext-delta))], [0.7 0.7 0.7], 'FaceAlpha',0.6, 'EdgeColor','none')
set(gca,'box','off')
xlabel("autism factor score");
ylabel("traitAnx factor score");
h = legend({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'});
title(h,'group')


%% additional plots for reference (not in paper)

%plot IQ per group
vals = recap_excl.ICAR_score;
mean_f = [mean(vals(grp==1)) mean(vals(grp==2)) mean(vals(grp==3)) mean(vals(grp==4)) mean(vals(grp==5))];
sem_f = [std(vals(grp==1))/sqrt(gsize(1)) std(vals(grp==2))/sqrt(gsize(2)) ...
    std(vals(grp==3))/sqrt(gsize(3)) std(vals(grp==4))/sqrt(gsize(4)) std(vals(grp==5))/sqrt(gsize(5))];
figure;
b = bar(1:5,mean_f,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
b.CData = clr; hold on
plotSpread(vals,'distributionIdx',grp,'distributionColors','k');
errorbar(1:5,mean_f,sem_f,'.k','LineWidth',1.5);
xticks(1:5)
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
ylabel('ICAR score')
xlabel('Group')
set(gca,'box','off')

%plot age per group
vals = recap_excl.age;
mean_f = [mean(vals(grp==1)) mean(vals(grp==2)) mean(vals(grp==3)) mean(vals(grp==4)) mean(vals(grp==5))];
sem_f = [std(vals(grp==1))/sqrt(gsize(1)) std(vals(grp==2))/sqrt(gsize(2)) ...
    std(vals(grp==3))/sqrt(gsize(3)) std(vals(grp==4))/sqrt(gsize(4)) std(vals(grp==5))/sqrt(gsize(5))];
figure;
b = bar(1:5,mean_f,0.7,'facecolor','flat','EdgeColor','k','LineWidth',1);
b.CData = clr; hold on
plotSpread(vals,'distributionIdx',grp,'distributionColors','k');
errorbar(1:5,mean_f,sem_f,'.k','LineWidth',1.5);
xticks(1:5)
xticklabels({'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'})
xtickangle(30)
ylim([18 65])
ylabel('Age')
xlabel('Group')
set(gca,'box','off')

%for comparison, plot mean questionnaire summary scores per group
g_list = {'Baseline','ExpLearn','ObsLearn','FixArb','DynArb'};
recap_Q = zscore(table2array(recap_excl(:,11:15)));
figure;
for g = 1:5
    vals = recap_Q(grp==g,:);
    subplot(2,3,g); hold on
    b = bar(1:5,mean(vals),0.7,'facecolor',clr(g,:),'EdgeColor','k','LineWidth',1);
    %plotSpread(vals,'distributionIdx',grp,'distributionColors','k');
    errorbar(1:5,mean(vals),std(vals)/sqrt(gsize(g)),'.k','LineWidth',1.5);
    xticks(1:5)
    xticklabels({'STAI State','STAI Trait','BDI','SRS','LSAS'})
    xtickangle(30)
    ylabel('factor score')
    title([g_list{g} ' group'])
end


%plot F3_autism and F7_traitAnxiety for Baseline and DynArb groups specifically
x = recap_excl.F3_autism;
y = recap_excl.F7_traitAnx;
i15 = recap_excl.group==1 | recap_excl.group==5;
figure;
scatterhist(x(i15),y(i15),'Group',recap_excl.group(i15),'Kernel','on',...
    'LineStyle','-','Color',clr([1 5],:),'Marker','.')
set(gca,'box','off')
xlabel("F3 score (~Autism traits)");
ylabel("F7 score (~Trait anxiety)");
legend({'Baseline group','DynArb group'})


%% correlation matrix between factor and parameters (not in paper)

corr_mat_f_p = table('Size',[8 10],'VariableTypes',...
    {'double','double','double','double','double','double','double','double','double','double'},...
    'VariableNames', {'factor','colorBias','handBias','stickyAct','actImit',...
    'EL_alpha','EL_magBoost','OL_alpha','wOLEL_fix','biasOLEL_dyn'});
for f=1:8
    corr_mat_f_p.factor(f) = f;
    for p=1:9
        factor_score = table2array(recap_excl(:,82+f));
        if p<=2
            parameter = abs(table2array(recap_excl(:,55+p))); %absolute value for color bias and hand bias
        else
            parameter = table2array(recap_excl(:,55+p));
        end
        [R,P] = corrcoef(parameter,factor_score);
        corr_mat_f_p(f,p+1) = array2table(R(1,2));
        p_values(f,p) = P(1,2);
    end
end
cdata = table2array(corr_mat_f_p(:,2:10));
xvalues = {'colorBias','handBias','stickyAct','actImit','ELalpha','ELmagBoost','OLalpha','wOLELfix','biasOLELdyn'};
yvalues = {'depression','socialAnx','autism','stateAnx','socialResp','groupAvoid','traitAnx','perfAnx'};
figure;
h = heatmap(xvalues,yvalues,cdata);
h.XLabel = 'Parameter';
h.YLabel = 'Factor';
h.Colormap = cool;

figure;
h = heatmap(xvalues,yvalues,p_values);
h.XLabel = 'Parameter';
h.YLabel = 'Factor';
h.Title = 'P-values';

corr(recap_excl.EL_alpha,recap_excl.Resp_ExpLearn)
corr(recap_excl.EL_magBoost,recap_excl.Resp_ExpLearn)
corr(recap_excl.OL_alpha,recap_excl.Resp_ObsLearn)
corr(recap_excl.EL_alpha,recap_excl.Resp_FixArb)

%same analysis with summary score to check that factor associations are stronger
corr_mat = nan(5,9);
p_values2 = nan(5,9);
for f=1:5
    for p=1:9
        q_score = table2array(recap_excl(:,10+f));
        if p<=2
            parameter = abs(table2array(recap_excl(:,55+p))); %absolute value for color bias and hand bias
        else
            parameter = table2array(recap_excl(:,55+p));
        end
        [R,P] = corrcoef(parameter,q_score);
        corr_mat(f,p) = R(1,2);
        p_values2(f,p) = P(1,2);
    end
end
xvalues = {'colorBias','handBias','stickyAct','actImit','ELalpha','ELmagBoost','OLalpha','wOLELfix','biasOLELdyn'};
yvalues = {'STAI State','STAI Trait','BDI','SRS','LSAS'};
figure;
h = heatmap(xvalues,yvalues,corr_mat);
h.XLabel = 'Parameter';
h.YLabel = 'Factor';
h.Colormap = cool;
h.ColorLimits = [-0.14 0.19];

figure;
h = heatmap(xvalues,yvalues,p_values2);
h.XLabel = 'Parameter';
h.YLabel = 'Factor';
h.Title = 'P-values';