%% 
% Author: Benoit Caldairou, PhD
% mail: benoit.caldairou@mcgill.ca
% This script performs a statistical study comparing the TLE and the
% controls. It can safely be ignored.


%% Surface-based statistical study -- Whole population

% ipsi-contra
groups_simplified = groups_training;
groups_simplified(tle_group_training) = {'TLE'};

z_columnar_volume  = [z_columvol_patients_ipsi_contra;z_columvol_controls];
z_t2_signal        = [z_t2_signal_patients_ipsi_contra;z_t2_signal_controls];
z_ratio_nuc_signal = [z_ratio_nuc_signal_patients_ipsi_contra;z_ratio_nuc_signal_controls];

% Build the model
Groups = term(cellstr(groups_simplified));
M = 1 + Groups;

% Fit the model
slm_columvol         = SurfStatLinMod(z_columnar_volume, M, template);
slm_t2_signal        = SurfStatLinMod(z_t2_signal, M, template);
slm_ratio_nuc_signal = SurfStatLinMod(z_ratio_nuc_signal, M, template);

% Contrast
pc_contrast = Groups.TLE - Groups.Controls;
cp_contrast = Groups.Controls - Groups.TLE;

% Get T maps
slm_columvol_up = SurfStatT(slm_columvol,pc_contrast);
slm_columvol_dn = SurfStatT(slm_columvol,cp_contrast);
slm_t2_signal = SurfStatT(slm_t2_signal,pc_contrast);
slm_ratio_nuc_signal = SurfStatT(slm_ratio_nuc_signal,pc_contrast);

figure(1), BladeSurfStatViewData(slm_columvol_up.t, template, 't map -- columnar volume -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_ColVol_IpsiContra_tMap.tif','-dtiffn','-r600');
figure(2), BladeSurfStatViewData(slm_t2_signal.t, template, 't map -- t2-signal -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_T2_IpsiContra_tMap.tif','-dtiffn','-r600');
figure(3), BladeSurfStatViewData(slm_ratio_nuc_signal.t, template, 't map -- ratio-nuc-signal -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_Ratio_IpsiContra_tMap.tif','-dtiffn','-r600');

% Run Model
[pval_columvol_dn, peak_columvol_dn, clus_columvol_dn] = SurfStatP(slm_columvol_dn, logical(ones(1,size(template.coord,2))), 0.05); 
[pval_t2_signal, peak_t2_signal, clus_t2_signal] = SurfStatP(slm_t2_signal, logical(ones(1,size(template.coord,2))), 0.05); 
[pval_ratio_nuc_signal, peak_ratio_nuc_signal, clus_ratio_nuc_signal] = SurfStatP(slm_ratio_nuc_signal, logical(ones(1,size(template.coord,2))), 0.05); 

% Uncorrected P-Values
% figure(7), BladeSurfStatViewData(1-tcdf(slm_columvol_dn.t, slm_columvol_dn.df), template, 'UnCorrected P-Values -- columnar volume -- ipsi/contra'); SurfStatColLim([0 0.05]);
% figure(8), BladeSurfStatViewData(1-tcdf(slm_t2_signal.t, slm_t2_signal.df), template, 'UnCorrected P-Values -- t2 signal -- ipsi/contra'); SurfStatColLim([0 0.05]);
% figure(9), BladeSurfStatViewData(1-tcdf(slm_ratio_nuc_signal.t, slm_ratio_nuc_signal.df), template, 'UnCorrected P-Values -- ratio NUC signal -- ipsi/contra'); SurfStatColLim([0 0.05]);

% Cohen Effect Size ColVol
D_CV = slm_columvol_up.ef;
D_THR_CV =D_CV .* ( pval_columvol_dn.C < 0.05);

figure(10),BladeSurfStatViewData(D_THR_CV, template, 'Cohens Effect Size Columnar Volume (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_ColVol_IpsiContra_EffectSize.tif','-dtiffn','-r600');

% Cohen Effect Size T2 Signal
D_intensity = slm_t2_signal.ef;
D_THR_intensity =D_intensity .* (pval_t2_signal.C < 0.05);

figure(11),BladeSurfStatViewData(D_THR_intensity, template, 'Cohens Effect Size T2 Signal (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_T2_IpsiContra_EffectSize.tif','-dtiffn','-r600');

% Cohen Effect Size Ratio NUC Signal
D_ratio_nuc = slm_ratio_nuc_signal.ef;
D_THR_ratio_nuc = D_ratio_nuc .* (pval_ratio_nuc_signal.C < 0.05);

figure(12),BladeSurfStatViewData(D_THR_ratio_nuc, template, 'Cohens Effect Size Ratio NUC Signal (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_Ratio_IpsiContra_EffectSize.tif','-dtiffn','-r600');

%% Surface-based statistical study -- Whole Population -- Assymetry

groups_simplified = groups_training;
groups_simplified(tle_group_training) = {'TLE'};

z_columnar_volume  = [z_columvol_ass_patients_ipsi_contra;z_columvol_ass_controls];
z_t2_signal        = [z_t2_signal_ass_patients_ipsi_contra;z_t2_signal_ass_controls];
z_ratio_nuc_signal = [z_ratio_nuc_signal_ass_patients_ipsi_contra;z_ratio_nuc_signal_ass_controls];

% Build the model
Groups = term(cellstr(groups_simplified));
M = 1 + Groups;

% Fit the model
slm_columvol       = SurfStatLinMod(z_columnar_volume, M, template_uni);
slm_t2_signal      = SurfStatLinMod(z_t2_signal, M, template_uni);
slm_ratio_nuc_signal   = SurfStatLinMod(z_ratio_nuc_signal, M, template_uni);

% Contrast
pc_contrast = Groups.TLE - Groups.Controls;
cp_contrast = Groups.Controls - Groups.TLE;

slm_columvol_up = SurfStatT(slm_columvol,pc_contrast);
slm_columvol_dn = SurfStatT(slm_columvol,cp_contrast);
slm_t2_signal = SurfStatT(slm_t2_signal,pc_contrast);
slm_ratio_nuc_signal = SurfStatT(slm_ratio_nuc_signal,pc_contrast);

% T Maps
figure(1), BladeSurfStatViewData([slm_columvol_up.t,zeros(1,21766)], template, 't map -- columnar volume -- ipsi/contra'); SurfStatColLim([-7 7]);
figure(2), BladeSurfStatViewData([slm_t2_signal.t,zeros(1,21766)], template, 't map -- t2-signal -- ipsi/contra'); SurfStatColLim([-7 7]);
figure(3), BladeSurfStatViewData([slm_ratio_nuc_signal.t,zeros(1,21766)], template, 't map -- ratio-nuc-signal -- ipsi/contra'); SurfStatColLim([-7 7]);

% Run Model
[pval_columvol_dn, peak_columvol_dn, clus_columvol_dn, clusid_columvol_dn] = SurfStatP(slm_columvol_dn, logical(ones(1,size(template_uni.coord,2))), 0.01); 
[pval_t2_signal, peak_t2_signal, clus_t2_signal] = SurfStatP(slm_t2_signal, logical(ones(1,size(template_uni.coord,2))), 0.01); 
[pval_ratio_nuc_signal, peak_ratio_nuc_signal, clus_ratio_nuc_signal] = SurfStatP(slm_ratio_nuc_signal, logical(ones(1,size(template_uni.coord,2))), 0.01); 

% Display uncorrected p-values
figure(4), BladeSurfStatViewData([1-tcdf(slm_columvol_dn.t, slm_columvol_dn.df),zeros(1,21766)], template, 'UnCorrected P-Values -- columnar volume -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(5), BladeSurfStatViewData([1-tcdf(slm_t2_signal.t, slm_t2_signal.df),zeros(1,21766)], template, 'UnCorrected P-Values -- t2 signal -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(6), BladeSurfStatViewData([1-tcdf(slm_ratio_nuc_signal.t, slm_ratio_nuc_signal.df),zeros(1,21766)], template, 'UnCorrected P-Values -- ratio nuc signal -- ipsi/contra'); SurfStatColLim([0 0.05]);

% Cohen Effect Size Volume
D_CV = slm_columvol_up.ef;
D_THR_CV =D_CV .* ( pval_columvol_dn.C < 0.05);

figure(7),BladeSurfStatViewData([D_THR_CV,zeros(1,21766)], template, 'Cohens Effect Size_columnar volume (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);

% Cohen Effect Size T2 Signal
D_intensity = slm_t2_signal.ef;
D_THR_intensity =D_intensity .* (pval_t2_signal.C < 0.05);

figure(8),BladeSurfStatViewData([D_THR_intensity,zeros(1,21766)], template, 'Cohens Effect Size T2 (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);

% Cohen Effect Size Ratio NUC Signal
D_ratio_nuc = slm_ratio_nuc_signal.ef;
D_THR_ratio_nuc =D_ratio_nuc .* (pval_ratio_nuc_signal.C < 0.05);

figure(9),BladeSurfStatViewData([D_THR_ratio_nuc,zeros(1,21766)], template, 'Cohens Effect Size Ratio NUC (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);


%% Surface-based statistical study -- HS population

groups_simplified = groups_training;
groups_simplified(tle_group_training) = {'TLE'};
groups_simplified(tle_group_training & (strcmp(hs_clinical_training,'no'))) = [];

% ipsi-contra
z_columnar_volume  = [z_columvol_patients_ipsi_contra;z_columvol_controls];
z_t2_signal        = [z_t2_signal_patients_ipsi_contra;z_t2_signal_controls];
z_ratio_nuc_signal = [z_ratio_nuc_signal_patients_ipsi_contra;z_ratio_nuc_signal_controls];

% Remove patients who are not HS
z_columnar_volume(tle_group_training & (strcmp(hs_clinical_training,'no')),:)  = [];
z_t2_signal(tle_group_training & (strcmp(hs_clinical_training,'no')),:)        = [];
z_ratio_nuc_signal(tle_group_training & (strcmp(hs_clinical_training,'no')),:) = [];

% Build the model
Groups = term(cellstr(groups_simplified));
M = 1 + Groups;

% Fit the model
slm_columvol         = SurfStatLinMod(z_columnar_volume, M, template);
slm_t2_signal        = SurfStatLinMod(z_t2_signal, M, template);
slm_ratio_nuc_signal = SurfStatLinMod(z_ratio_nuc_signal, M, template);

% Contrast
pc_contrast = Groups.TLE - Groups.Controls;
cp_contrast = Groups.Controls - Groups.TLE;

% Get T maps
slm_columvol_up = SurfStatT(slm_columvol,pc_contrast);
slm_columvol_dn = SurfStatT(slm_columvol,cp_contrast);
slm_t2_signal = SurfStatT(slm_t2_signal,pc_contrast);
slm_ratio_nuc_signal = SurfStatT(slm_ratio_nuc_signal,pc_contrast);

figure(1), BladeSurfStatViewData(slm_columvol_up.t, template, 't map -- columnar volume -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_HS_ColVol_IpsiContra_tMap.tif','-dtiffn','-r600');
figure(2), BladeSurfStatViewData(slm_t2_signal.t, template, 't map -- t2-signal -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_HS_T2_IpsiContra_tMap.tif','-dtiffn','-r600');
figure(3), BladeSurfStatViewData(slm_ratio_nuc_signal.t, template, 't map -- ratio-nuc-signal -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_HS_Ratio_IpsiContra_tMap.tif','-dtiffn','-r600');

% Run Model
[pval_columvol_dn, peak_columvol_dn, clus_columvol_dn] = SurfStatP(slm_columvol_dn, logical(ones(1,size(template.coord,2))), 0.05); 
[pval_t2_signal, peak_t2_signal, clus_t2_signal] = SurfStatP(slm_t2_signal, logical(ones(1,size(template.coord,2))), 0.05); 
[pval_ratio_nuc_signal, peak_ratio_nuc_signal, clus_ratio_nuc_signal] = SurfStatP(slm_ratio_nuc_signal, logical(ones(1,size(template.coord,2))), 0.05); 

% Uncorrected P-Values
figure(4), BladeSurfStatViewData(1-tcdf(slm_columvol_dn.t, slm_columvol_dn.df), template, 'UnCorrected P-Values -- columnar volume -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(5), BladeSurfStatViewData(1-tcdf(slm_t2_signal.t, slm_t2_signal.df), template, 'UnCorrected P-Values -- t2 signal -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(6), BladeSurfStatViewData(1-tcdf(slm_ratio_nuc_signal.t, slm_ratio_nuc_signal.df), template, 'UnCorrected P-Values -- ratio NUC signal -- ipsi/contra'); SurfStatColLim([0 0.05]);

% Cohen Effect Size ColVol
D_CV = slm_columvol_up.ef;
D_THR_CV =D_CV .* ( pval_columvol_dn.C < 0.05);

figure(7),BladeSurfStatViewData(D_THR_CV, template, 'Cohens Effect Size Columnar Volume (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_HS_ColVol_IpsiContra_EffectSize.tif','-dtiffn','-r600');

% Cohen Effect Size T2 Signal
D_intensity = slm_t2_signal.ef;
D_THR_intensity =D_intensity .* (pval_t2_signal.C < 0.05);

figure(8),BladeSurfStatViewData(D_THR_intensity, template, 'Cohens Effect Size T2 Signal (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_HS_T2_IpsiContra_EffectSize.tif','-dtiffn','-r600');

% Cohen Effect Size Ratio NUC Signal
D_ratio_nuc = slm_ratio_nuc_signal.ef;
D_THR_ratio_nuc = D_ratio_nuc .* (pval_ratio_nuc_signal.C < 0.05);

figure(9),BladeSurfStatViewData(D_THR_ratio_nuc, template, 'Cohens Effect Size Ratio NUC Signal (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_HS_Ratio_IpsiContra_EffectSize.tif','-dtiffn','-r600');

%% Surface-based statistical study -- HS Population -- Assymetry

groups_simplified = groups_training;
groups_simplified(tle_group_training) = {'TLE'};
groups_simplified(tle_group_training & (strcmp(hs_clinical_training,'no'))) = [];

% ipsi-contra
z_columnar_volume  = [z_columvol_ass_patients_ipsi_contra;z_columvol_ass_controls];
z_t2_signal        = [z_t2_signal_ass_patients_ipsi_contra;z_t2_signal_ass_controls];
z_ratio_nuc_signal = [z_ratio_nuc_signal_ass_patients_ipsi_contra;z_ratio_nuc_signal_ass_controls];

% Remove patients who are not HS
z_columnar_volume(tle_group_training & (strcmp(hs_clinical_training,'no')),:) = [];
z_t2_signal(tle_group_training & (strcmp(hs_clinical_training,'no')),:) = [];
z_ratio_nuc_signal(tle_group_training & (strcmp(hs_clinical_training,'no')),:) = [];

% Build the model
Groups = term(cellstr(groups_simplified));
M = 1 + Groups;

% Fit the model
slm_columvol       = SurfStatLinMod(z_columnar_volume, M, template_uni);
slm_t2_signal      = SurfStatLinMod(z_t2_signal, M, template_uni);
slm_ratio_nuc_signal   = SurfStatLinMod(z_ratio_nuc_signal, M, template_uni);

% Contrast
pc_contrast = Groups.TLE - Groups.Controls;
cp_contrast = Groups.Controls - Groups.TLE;
slm_columvol_up = SurfStatT(slm_columvol,pc_contrast);
slm_columvol_dn = SurfStatT(slm_columvol,cp_contrast);
slm_t2_signal = SurfStatT(slm_t2_signal,pc_contrast);
slm_ratio_nuc_signal = SurfStatT(slm_ratio_nuc_signal,pc_contrast);

% T Maps
figure(1), BladeSurfStatViewData([slm_columvol_up.t,zeros(1,21766)], template, 't map -- columnar volume -- ipsi/contra'); SurfStatColLim([-7 7]);
figure(2), BladeSurfStatViewData([slm_t2_signal.t,zeros(1,21766)], template, 't map -- t2-signal -- ipsi/contra'); SurfStatColLim([-7 7]);
figure(3), BladeSurfStatViewData([slm_ratio_nuc_signal.t,zeros(1,21766)], template, 't map -- ratio-nuc-signal -- ipsi/contra'); SurfStatColLim([-7 7]);

% Run Model
[pval_columvol_dn, peak_columvol_dn, clus_columvol_dn, clusid_columvol_dn] = SurfStatP(slm_columvol_dn, logical(ones(1,size(template_uni.coord,2))), 0.01); 
[pval_t2_signal, peak_t2_signal, clus_t2_signal] = SurfStatP(slm_t2_signal, logical(ones(1,size(template_uni.coord,2))), 0.01); 
[pval_ratio_nuc_signal, peak_ratio_nuc_signal, clus_ratio_nuc_signal] = SurfStatP(slm_ratio_nuc_signal, logical(ones(1,size(template_uni.coord,2))), 0.01); 

% Display uncorrected p-values
figure(4), BladeSurfStatViewData([1-tcdf(slm_columvol_dn.t, slm_columvol_dn.df),zeros(1,21766)], template, 'UnCorrected P-Values -- columnar volume -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(5), BladeSurfStatViewData([1-tcdf(slm_t2_signal.t, slm_t2_signal.df),zeros(1,21766)], template, 'UnCorrected P-Values -- t2 signal -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(6), BladeSurfStatViewData([1-tcdf(slm_ratio_nuc_signal.t, slm_ratio_nuc_signal.df),zeros(1,21766)], template, 'UnCorrected P-Values -- ratio nuc signal -- ipsi/contra'); SurfStatColLim([0 0.05]);

% Cohen Effect Size Volume
D_CV = slm_columvol_up.ef;
D_THR_CV =D_CV .* ( pval_columvol_dn.C < 0.05);

figure(7),BladeSurfStatViewData([D_THR_CV,zeros(1,21766)], template, 'Cohens Effect Size_columnar volume (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);

% Cohen Effect Size T2 Signal
D_intensity = slm_t2_signal.ef;
D_THR_intensity =D_intensity .* (pval_t2_signal.C < 0.05);

figure(8),BladeSurfStatViewData([D_THR_intensity,zeros(1,21766)], template, 'Cohens Effect Size T2 (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);

% Cohen Effect Size Ratio NUC Signal
D_ratio_nuc = slm_ratio_nuc_signal.ef;
D_THR_ratio_nuc =D_ratio_nuc .* (pval_ratio_nuc_signal.C < 0.05);

figure(9),BladeSurfStatViewData([D_THR_ratio_nuc,zeros(1,21766)], template, 'Cohens Effect Size Ratio NUC (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);

%% Surface-based statistical study -- nonHS population

groups_simplified = groups_training;
groups_simplified(tle_group_training) = {'TLE'};
groups_simplified(tle_group_training & (strcmp(hs_clinical_training,'yes'))) = [];

% ipsi-contra
z_columnar_volume  = [z_columvol_patients_ipsi_contra;z_columvol_controls];
z_t2_signal        = [z_t2_signal_patients_ipsi_contra;z_t2_signal_controls];
z_ratio_nuc_signal = [z_ratio_nuc_signal_patients_ipsi_contra;z_ratio_nuc_signal_controls];

% Remove patients who are HS
z_columnar_volume(tle_group & (strcmp(hs_clinical_training,'yes')),:)  = [];
z_t2_signal(tle_group & (strcmp(hs_clinical_training,'yes')),:)        = [];
z_ratio_nuc_signal(tle_group & (strcmp(hs_clinical_training,'yes')),:) = [];

% Build the model
Groups = term(cellstr(groups_simplified));
M = 1 + Groups;

% Fit the model
slm_columvol         = SurfStatLinMod(z_columnar_volume, M, template);
slm_t2_signal        = SurfStatLinMod(z_t2_signal, M, template);
slm_ratio_nuc_signal = SurfStatLinMod(z_ratio_nuc_signal, M, template);

% Contrast
pc_contrast = Groups.TLE - Groups.Controls;
cp_contrast = Groups.Controls - Groups.TLE;

% Get T maps
slm_columvol_up = SurfStatT(slm_columvol,pc_contrast);
slm_columvol_dn = SurfStatT(slm_columvol,cp_contrast);
slm_t2_signal = SurfStatT(slm_t2_signal,pc_contrast);
slm_ratio_nuc_signal = SurfStatT(slm_ratio_nuc_signal,pc_contrast);

figure(1), BladeSurfStatViewData(slm_columvol_up.t, template, 't map -- columnar volume -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_nonHS_ColVol_IpsiContra_tMap.tif','-dtiffn','-r600');
figure(2), BladeSurfStatViewData(slm_t2_signal.t, template, 't map -- t2-signal -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_nonHS_T2_IpsiContra_tMap.tif','-dtiffn','-r600');
figure(3), BladeSurfStatViewData(slm_ratio_nuc_signal.t, template, 't map -- ratio-nuc-signal -- ipsi/contra'); SurfStatColLim([-7 7]);print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_nonHS_Ratio_IpsiContra_tMap.tif','-dtiffn','-r600');

% Run Model
[pval_columvol_dn, peak_columvol_dn, clus_columvol_dn] = SurfStatP(slm_columvol_dn, logical(ones(1,size(template.coord,2))), 0.05); 
[pval_t2_signal, peak_t2_signal, clus_t2_signal] = SurfStatP(slm_t2_signal, logical(ones(1,size(template.coord,2))), 0.05); 
[pval_ratio_nuc_signal, peak_ratio_nuc_signal, clus_ratio_nuc_signal] = SurfStatP(slm_ratio_nuc_signal, logical(ones(1,size(template.coord,2))), 0.05); 

% Uncorrected P-Values
figure(4), BladeSurfStatViewData(1-tcdf(slm_columvol_dn.t, slm_columvol_dn.df), template, 'UnCorrected P-Values -- columnar volume -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(5), BladeSurfStatViewData(1-tcdf(slm_t2_signal.t, slm_t2_signal.df), template, 'UnCorrected P-Values -- t2 signal -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(6), BladeSurfStatViewData(1-tcdf(slm_ratio_nuc_signal.t, slm_ratio_nuc_signal.df), template, 'UnCorrected P-Values -- ratio NUC signal -- ipsi/contra'); SurfStatColLim([0 0.05]);

% Cohen Effect Size ColVol
D_CV = slm_columvol_up.ef;
D_THR_CV =D_CV .* ( pval_columvol_dn.C < 0.05);

figure(7),BladeSurfStatViewData(D_THR_CV, template, 'Cohens Effect Size Columnar Volume (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_nonHS_ColVol_IpsiContra_EffectSize.tif','-dtiffn','-r600');

% Cohen Effect Size T2 Signal
D_intensity = slm_t2_signal.ef;
D_THR_intensity =D_intensity .* (pval_t2_signal.C < 0.05);

figure(8),BladeSurfStatViewData(D_THR_intensity, template, 'Cohens Effect Size T2 Signal (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_nonHS_T2_IpsiContra_EffectSize.tif','-dtiffn','-r600');

% Cohen Effect Size Ratio NUC Signal
D_ratio_nuc = slm_ratio_nuc_signal.ef;
D_THR_ratio_nuc = D_ratio_nuc .* (pval_ratio_nuc_signal.C < 0.05);

figure(9),BladeSurfStatViewData(D_THR_ratio_nuc, template, 'Cohens Effect Size Ratio NUC Signal (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
print(gcf,'ForPPT_SeparateHistopatho/All_nonHisto_nonHS_Ratio_IpsiContra_EffectSize.tif','-dtiffn','-r600');

%% Quick surface-based statistical study -- nonHS Population -- Assymetry

groups_simplified = groups_training;
groups_simplified(tle_group_training) = {'TLE'};
groups_simplified(tle_group_training & (strcmp(hs_clinical,'yes'))) = [];

% ipsi-contra
z_columnar_volume  = [z_columvol_ass_patients_ipsi_contra;z_columvol_ass_controls];
z_t2_signal        = [z_t2_signal_ass_patients_ipsi_contra;z_t2_signal_ass_controls];
z_ratio_nuc_signal = [z_ratio_nuc_signal_ass_patients_ipsi_contra;z_ratio_nuc_signal_ass_controls];

% Remove patients who are HS
z_columnar_volume(tle_group_training & (strcmp(hs_clinical_training,'yes')),:)  = [];
z_t2_signal(tle_group_training & (strcmp(hs_clinical_training,'yes')),:)        = [];
z_ratio_nuc_signal(tle_group_training & (strcmp(hs_clinical_training,'yes')),:) = [];

% Build the model
Groups = term(cellstr(groups_simplified));
M = 1 + Groups;

% Fit the model
slm_columvol       = SurfStatLinMod(z_columnar_volume, M, template_uni);
slm_t2_signal      = SurfStatLinMod(z_t2_signal, M, template_uni);
slm_ratio_nuc_signal   = SurfStatLinMod(z_ratio_nuc_signal, M, template_uni);

% Contrast
pc_contrast = Groups.TLE - Groups.Controls;
cp_contrast = Groups.Controls - Groups.TLE;
slm_columvol_up = SurfStatT(slm_columvol,pc_contrast);
slm_columvol_dn = SurfStatT(slm_columvol,cp_contrast);
slm_t2_signal = SurfStatT(slm_t2_signal,pc_contrast);
slm_ratio_nuc_signal = SurfStatT(slm_ratio_nuc_signal,pc_contrast);

% T Maps
figure(1), BladeSurfStatViewData([slm_columvol_up.t,zeros(1,21766)], template, 't map -- columnar volume -- ipsi/contra'); SurfStatColLim([-7 7]);
figure(2), BladeSurfStatViewData([slm_t2_signal.t,zeros(1,21766)], template, 't map -- t2-signal -- ipsi/contra'); SurfStatColLim([-7 7]);
figure(3), BladeSurfStatViewData([slm_ratio_nuc_signal.t,zeros(1,21766)], template, 't map -- ratio-nuc-signal -- ipsi/contra'); SurfStatColLim([-7 7]);

% Run Model
[pval_columvol_dn, peak_columvol_dn, clus_columvol_dn, clusid_columvol_dn] = SurfStatP(slm_columvol_dn, logical(ones(1,size(template_uni.coord,2))), 0.01); 
[pval_t2_signal, peak_t2_signal, clus_t2_signal] = SurfStatP(slm_t2_signal, logical(ones(1,size(template_uni.coord,2))), 0.01); 
[pval_ratio_nuc_signal, peak_ratio_nuc_signal, clus_ratio_nuc_signal] = SurfStatP(slm_ratio_nuc_signal, logical(ones(1,size(template_uni.coord,2))), 0.01); 

% Display uncorrected p-values
figure(4), BladeSurfStatViewData([1-tcdf(slm_columvol_dn.t, slm_columvol_dn.df),zeros(1,21766)], template, 'UnCorrected P-Values -- columnar volume -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(5), BladeSurfStatViewData([1-tcdf(slm_t2_signal.t, slm_t2_signal.df),zeros(1,21766)], template, 'UnCorrected P-Values -- t2 signal -- ipsi/contra'); SurfStatColLim([0 0.05]);
figure(6), BladeSurfStatViewData([1-tcdf(slm_ratio_nuc_signal.t, slm_ratio_nuc_signal.df),zeros(1,21766)], template, 'UnCorrected P-Values -- ratio nuc signal -- ipsi/contra'); SurfStatColLim([0 0.05]);

% Cohen Effect Size Volume
D_CV = slm_columvol_up.ef;
D_THR_CV =D_CV .* ( pval_columvol_dn.C < 0.05);

figure(7),BladeSurfStatViewData([D_THR_CV,zeros(1,21766)], template, 'Cohens Effect Size_columnar volume (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);

% Cohen Effect Size T2 Signal
D_intensity = slm_t2_signal.ef;
D_THR_intensity =D_intensity .* (pval_t2_signal.C < 0.05);

figure(8),BladeSurfStatViewData([D_THR_intensity,zeros(1,21766)], template, 'Cohens Effect Size T2 (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);

% Cohen Effect Size Ratio NUC Signal
D_ratio_nuc = slm_ratio_nuc_signal.ef;
D_THR_ratio_nuc =D_ratio_nuc .* (pval_ratio_nuc_signal.C < 0.05);

figure(9),BladeSurfStatViewData([D_THR_ratio_nuc,zeros(1,21766)], template, 'Cohens Effect Size Ratio NUC (P_F_W_E < 0.05)');
SurfStatColLim([-1 1]);
colormap(cohen_colormap);
