% Statistical analysis of all cohort data

load('t1_RAREVTR_AllData.mat');
genotypes = t1RAREVTRAllData.Genotype;
timepoints = t1RAREVTRAllData.Timepoint;
infected = t1RAREVTRAllData.Infected;
emd_s = t1RAREVTRAllData.EMDstriatum;
emd_t = t1RAREVTRAllData.EMDthalamus;
emd_p = t1RAREVTRAllData.EMDpssc;
emd_d = t1RAREVTRAllData.EMDdg;

% timepoint 1
genotypes_1 = genotypes(timepoints == 1);
infected_1 = categorical(infected(timepoints == 1));
emd_s1 = emd_s(timepoints == 1);
emd_t1 = emd_t(timepoints == 1);
emd_p1 = emd_p(timepoints == 1);
emd_d1 = emd_d(timepoints == 1);

wtni = emd_d1(genotypes_1 == 'WT' & infected_1 == '0');
wtif = emd_d1(genotypes_1 == 'WT' & infected_1 == '1');
koni = emd_d1(genotypes_1 == 'KO' & infected_1 == '0');
koif = emd_d1(genotypes_1 == 'KO' & infected_1 == '1');

genotype_anova = {'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; ...
                  'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF';...
                  'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; ...
                  'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF';};
close all; clc; 

figure(1); set(gcf, 'Color', 'w')
violinplot([wtni; wtif; koni; koif], genotype_anova, ...
                       'GroupOrder', {'WT-NI'; 'WT-IF'; 'KO-NI'; 'KO-IF'});              
ylabel('EMD (msec)');
title('Timepoint 1 - Dentate Gyrus');
              
% GROUPWISE ANOVA
figure()
[~, ~, stats] = anovan([wtni; wtif; koni; koif], ...
                        {genotype_anova});
anova_result = multcompare(stats)      

% Interaction ANOVA
[~,~,stats] = anovan(emd_d1,{genotypes_1, infected_1},'model','interaction',...
    'varnames',{'Genotype','Infection'});


%% timepoint 2
genotypes_2 = genotypes(timepoints == 2);
infected_2 = categorical(infected(timepoints == 2));
emd_s2 = emd_s(timepoints == 2);
emd_t2 = emd_t(timepoints == 2);
emd_p2 = emd_p(timepoints == 2);
emd_d2 = emd_d(timepoints == 2);


wtni = emd_d2(genotypes_2 == 'WT' & infected_2 == '0');
wtif = emd_d2(genotypes_2 == 'WT' & infected_2 == '1');
koni = emd_d2(genotypes_2 == 'KO' & infected_2 == '0');
koif = emd_d2(genotypes_2 == 'KO' & infected_2 == '1');

genotype_anova = {'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; 'WT-NI'; ...
                  'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; 'WT-IF'; ...
                  'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; 'KO-NI'; ...
                  'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF'; 'KO-IF';};
close all; clc; 

figure(1); set(gcf, 'Color', 'w')
violinplot([wtni; wtif; koni; koif], genotype_anova, ...
                       'GroupOrder', {'WT-NI'; 'WT-IF'; 'KO-NI'; 'KO-IF'});              
ylabel('EMD (msec)');
title('Timepoint 2 - Dentate Gyrus');
              
% GROUPWISE ANOVA
figure()
[~, ~, stats] = anovan([wtni; wtif; koni; koif], ...
                        {genotype_anova});
anova_result = multcompare(stats)      

% Interaction ANOVA
[~,~,stats] = anovan(emd_d2,{genotypes_2, infected_2},'model','interaction',...
    'varnames',{'Genotype','Infection'});


