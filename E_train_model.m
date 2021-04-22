%% Author: Benoit Caldairou, PhD
%% Mail: benoit.caldairou@mcgill.ca

%% Set-up basic variables
type = 'diagquadratic';

% ------------------------- ColVol ------------------------------ %
% t-test to get t-maps beween ipsi and contra for training set
[~,~,~,stats_nohs] = ttest2(z_columvol_patients_ipsi_contra(:,1:size(z_columvol_patients,2)/2),z_columvol_patients_ipsi_contra(:,size(z_columvol_patients,2)/2+1:end), 0.025, 'both','unequal');
tmap = abs(stats_nohs.tstat);

% Looks for optimal t in training set
t_optimal = optimal_thres_lookup(tmap,z_columvol_ass_patients,lateralization_bin_training);

% Build data for this specific t
roi_colvol_final = tmap >= t_optimal;
database = mean(z_columvol_ass_patients(:,roi_colvol_final),2);

% Build the model and classify
lda_colvol_final = fitcdiscr(database,lateralization_bin_training,'DiscrimType',type);

% ------------------------ T2 SIGNAL ---------------------------- %
% t-test to get t-maps beween ipsi and contra for training set
[~,~,~,stats_nohs] = ttest2(z_t2_signal_patients_ipsi_contra(:,1:size(z_t2_signal_patients,2)/2),z_t2_signal_patients_ipsi_contra(:,size(z_t2_signal_patients,2)/2+1:end), 0.025, 'both','unequal');
tmap = abs(stats_nohs.tstat);

% Looks for optimal t in training set
t_optimal = optimal_thres_lookup(tmap,z_t2_signal_ass_patients,lateralization_bin_training);

% Build data for this specific t
roi_t2_final = tmap >= t_optimal;
database = mean(z_t2_signal_ass_patients(:,roi_t2_final),2);

% Build the model
lda_t2_final = fitcdiscr(database,lateralization_bin_training,'DiscrimType',type);

% -------------------- RATIO NUC SIGNAL ------------------------- %
% t-test to get t-maps beween ipsi and contra for training set
[~,~,~,stats_nohs] = ttest2(z_ratio_nuc_signal_patients_ipsi_contra(:,1:size(z_ratio_nuc_signal_patients,2)/2),z_ratio_nuc_signal_patients_ipsi_contra(:,size(z_ratio_nuc_signal_patients,2)/2+1:end), 0.025, 'both','unequal');
tmap = abs(stats_nohs.tstat);

% Looks for optimal t in training set
t_optimal = optimal_thres_lookup(tmap,z_ratio_nuc_signal_ass_patients,lateralization_bin_training);

% Build data for this specific t
roi_ratio_final = tmap >= t_optimal;
database = mean(z_ratio_nuc_signal_ass_patients(:,roi_ratio_final),2);

% Build the model
lda_ratio_final = fitcdiscr(database,lateralization_bin_training,'DiscrimType',type);

% ------------- MULTIVARIATE T2 Ratio ----------------------- %
clear forT2test_z_ipsi_contra_columnarvolume_training forT2test_z_ipsi_contra_ratiosignal_training forT2test_z_ipsi_contra_all_training;

% t-test to get t-maps beween ipsi and contra for training set
forT2test_z_ipsi_contra_t2_signal_training = [z_t2_signal_patients_ipsi_contra(:,1:size(z_columvol_patients,2)/2);z_t2_signal_patients_ipsi_contra(:,size(z_t2_signal_patients,2)/2+1:end)];
forT2test_z_ipsi_contra_ratiosignalnuc_training = [z_ratio_nuc_signal_patients_ipsi_contra(:,1:size(z_ratio_nuc_signal_patients,2)/2);z_ratio_nuc_signal_patients_ipsi_contra(:,size(z_ratio_nuc_signal_patients,2)/2+1:end)];
forT2test_z_ipsi_contra_all_training = cat(3,forT2test_z_ipsi_contra_t2_signal_training,forT2test_z_ipsi_contra_ratiosignalnuc_training);

groups_ipsi_contra = {};
groups_ipsi_contra(1:size(forT2test_z_ipsi_contra_t2_signal_training,1)/2) = {'IPSI'};
groups_ipsi_contra(size(forT2test_z_ipsi_contra_t2_signal_training,1)/2+1:size(forT2test_z_ipsi_contra_t2_signal_training,1)) = {'CONTRA'};
GROUPS_IPSI_CONTRA = term(cellstr(groups_ipsi_contra));

% Put a very simple model
M = 1 + GROUPS_IPSI_CONTRA;
% Hotelling T2 to get T-maps beween ipsi and contra
slm_all = SurfStatLinMod(forT2test_z_ipsi_contra_all_training, M, template_uni);
slm_all = SurfStatT(slm_all,GROUPS_IPSI_CONTRA.IPSI-GROUPS_IPSI_CONTRA.CONTRA);% A few more variables

tmap = abs(slm_all.t);

clear training_data database;

% Looks for optimal t in training set
training_data(:,:,1)  = z_t2_signal_ass_patients;
training_data(:,:,2)  = z_ratio_nuc_signal_ass_patients;
t_optimal = optimal_thres_lookup(tmap,training_data,lateralization_bin_training);

% Build data for this specific t
roi_t2_ratio_final = tmap >= t_optimal;
database(:,1) = mean(z_t2_signal_ass_patients(:,roi_t2_ratio_final),2);
database(:,2) = mean(z_ratio_nuc_signal_ass_patients(:,roi_t2_ratio_final),2);

% Build the model
lda_t2_ratio_final = fitcdiscr(database,lateralization_bin_training,'DiscrimType','linear');

% Save the rois and the models in a .mat file
save(['./','final_models'],'roi_colvol_final','lda_colvol_final','roi_t2_final','lda_t2_final','roi_ratio_final','lda_ratio_final','roi_t2_ratio_final','lda_t2_ratio_final');
