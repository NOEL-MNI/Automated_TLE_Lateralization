%% Author: Benoit Caldairou, PhD
%% Mail: benoit.caldairou@mcgill.ca

%% Set-up basic variables
feature_dir_individual='/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling';
image_dir_individual='';

%% Load pretrained models, roi and controls feature average and std for z-score
load controls_stats.mat;
load final_models.mat;

%% Read the data (Columnar Volume)

ca_columvol_left_file_individual   = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_CA_ColVol.txt');
ca_columvol_right_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_R_CA_ColVol.txt');
dg_columvol_left_file_individual   = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_DG_ColVol.txt');
dg_columvol_right_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_R_DG_ColVol.txt');
sub_columvol_left_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_SUB_ColVol.txt');
sub_columvol_right_file_individual = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individua,'_R_SUB_ColVol.txt');

ca_columvol_left_individual   = zeros(length(ids_individual),10242);
ca_columvol_right_individual  = zeros(length(ids_individual),10242);
dg_columvol_left_individual   = zeros(length(ids_individual),5762);
dg_columvol_right_individual  = zeros(length(ids_individual),5762);
sub_columvol_left_individual  = zeros(length(ids_individual),5762);
sub_columvol_right_individual = zeros(length(ids_individual),5762);

for i=1:length(ids_individual)

    ca_columvol_left_individual(i,:) = SurfStatReadData(ca_columvol_left_file_individual{i});
    ca_columvol_right_individual(i,:) = SurfStatReadData(ca_columvol_right_file_individual{i});
    dg_columvol_left_individual(i,:) = SurfStatReadData(dg_columvol_left_file_individual{i});
    dg_columvol_right_individual(i,:) = SurfStatReadData(dg_columvol_right_file_individual{i});
    sub_columvol_left_individual(i,:) = SurfStatReadData(sub_columvol_left_file_individual{i});
    sub_columvol_right_individual(i,:) = SurfStatReadData(sub_columvol_right_file_individual{i});

end

%% Read and normalize the data (T2 signal)

ca_t2_left_file_individual   = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_CA_nnt2.txt');
ca_t2_right_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_R_CA_nnt2.txt');
dg_t2_left_file_individual   = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_DG_nnt2.txt');
dg_t2_right_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_R_DG_nnt2.txt');
sub_t2_left_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_SUB_nnt2.txt');
sub_t2_right_file_individual = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_R_SUB_nnt2.txt');
norm_t2_file_individual      = strcat(image_dir_individual,'/',prefix_individual,'_',ids_individual,'_T2norm.mat');

ca_t2_left_individual   = zeros(length(ids_individual),10242);
ca_t2_right_individual  = zeros(length(ids_individual),10242);
dg_t2_left_individual   = zeros(length(ids_individual),5762);
dg_t2_right_individual  = zeros(length(ids_individual),5762);
sub_t2_left_individual  = zeros(length(ids_individual),5762);
sub_t2_right_individual = zeros(length(ids_individual),5762);

for i=1:length(ids_individual)

    ca_t2_left_individual(i,:)   = SurfStatReadData(ca_t2_left_file_individual{i});
    ca_t2_right_individual(i,:)  = SurfStatReadData(ca_t2_right_file_individual{i});
    dg_t2_left_individual(i,:)   = SurfStatReadData(dg_t2_left_file_individual{i});
    dg_t2_right_individual(i,:)  = SurfStatReadData(dg_t2_right_file_individual{i});
    sub_t2_left_individual(i,:)  = SurfStatReadData(sub_t2_left_file_individual{i});
    sub_t2_right_individual(i,:) = SurfStatReadData(sub_t2_right_file_individual{i});
    
    tmp = load(norm_t2_file_individual{i});
    
    ca_t2_left_individual(i,:)   = ca_t2_left_individual(i,:)./tmp.Int_ref;
    ca_t2_right_individual(i,:)  = ca_t2_right_individual(i,:)./tmp.Int_ref;
    dg_t2_left_individual(i,:)   = dg_t2_left_individual(i,:)./tmp.Int_ref;
    dg_t2_right_individual(i,:)  = dg_t2_right_individual(i,:)./tmp.Int_ref;
    sub_t2_left_individual(i,:)  = sub_t2_left_individual(i,:)./tmp.Int_ref;
    sub_t2_right_individual(i,:) = sub_t2_right_individual(i,:)./tmp.Int_ref;

end

%% Read and normalize the ratio NUC data
ca_ratio_nuc_left_file_individual   = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_CA_t2wt1wratio.txt');
ca_ratio_nuc_right_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_R_CA_t2wt1wratio.txt');
dg_ratio_nuc_left_file_individual   = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_DG_t2wt1wratio.txt');
dg_ratio_nuc_right_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_R_DG_t2wt1wratio.txt');
sub_ratio_nuc_left_file_individual  = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_L_SUB_t2wt1wratio.txt');
sub_ratio_nuc_right_file_individual = strcat(feature_dir_individual,'/',prefix_individual,'_',ids_individual,'_R_SUB_t2wt1wratio.txt');

ca_ratio_nuc_left_individual   = zeros(length(ids_individual),10242);
ca_ratio_nuc_right_individual  = zeros(length(ids_individual),10242);
dg_ratio_nuc_left_individual   = zeros(length(ids_individual),5762);
dg_ratio_nuc_right_individual  = zeros(length(ids_individual),5762);
sub_ratio_nuc_left_individual  = zeros(length(ids_individual),5762);
sub_ratio_nuc_right_individual = zeros(length(ids_individual),5762);

for i=1:length(ids_individual)

    ca_ratio_nuc_left_individual(i,:)   = SurfStatReadData(ca_ratio_nuc_left_file_individual{i});
    ca_ratio_nuc_right_individual(i,:)  = SurfStatReadData(ca_ratio_nuc_right_file_individual{i});
    dg_ratio_nuc_left_individual(i,:)   = SurfStatReadData(dg_ratio_nuc_left_file_individual{i});
    dg_ratio_nuc_right_individual(i,:)  = SurfStatReadData(dg_ratio_nuc_right_file_individual{i});
    sub_ratio_nuc_left_individual(i,:)  = SurfStatReadData(sub_ratio_nuc_left_file_individual{i});
    sub_ratio_nuc_right_individual(i,:) = SurfStatReadData(sub_ratio_nuc_right_file_individual{i});

end

%% z-score with respect to controls

% Global local volumes
columvol_individual               = [ca_columvol_left_individual,sub_columvol_left_individual,dg_columvol_left_individual,ca_columvol_right_individual,sub_columvol_right_individual,dg_columvol_right_individual];
columvol_individual               = SurfStatSmooth(columvol_individual, template, 3);
z_columvol_patients_individual    = (columvol_individual - repmat(mu_columvol_controls,[size(columvol_individual,1),1]))./repmat(std_columvol_controls,[size(columvol_individual,1),1]);

% T2 signal
t2_signal_individual              = [ca_t2_left_individual,sub_t2_left_individual,dg_t2_left_individual,ca_t2_right_individual,sub_t2_right_individual,dg_t2_right_individual];
t2_signal_individual              = SurfStatSmooth(t2_signal_individual, template, 3);
z_t2_signal_patients_individual   = (t2_signal_individual - repmat(mu_t2_signal_controls,[size(t2_signal_individual,1),1]))./repmat(std_t2_signal_controls,[size(t2_signal_individual,1),1]);

% Ratio NUC Signal
ratio_nuc_signal_individual               = [ca_ratio_nuc_left_individual,sub_ratio_nuc_left_individual,dg_ratio_nuc_left_individual,ca_ratio_nuc_right_individual,sub_ratio_nuc_right_individual,dg_ratio_nuc_right_individual];
ratio_nuc_signal_individual               = SurfStatSmooth(ratio_nuc_signal_individual, template, 3);
z_ratio_nuc_signal_patients_individual    = (ratio_nuc_signal_individual - repmat(mu_ratio_nuc_signal_controls,[size(ratio_nuc_signal_individual,1),1]))./repmat(std_ratio_nuc_signal_controls,[size(ratio_nuc_signal_individual,1),1]);

% Assymetry columnar volumes
columvol_assymetry_individual         = 2*(columvol_individual(:,1:size(columvol_individual,2)/2) - columvol_individual(:,size(columvol_individual,2)/2+1:end))./(columvol_individual(:,1:size(columvol_individual,2)/2) + columvol_individual(:,size(columvol_individual,2)/2+1:end));
z_columvol_ass_patients_individual    = (columvol_assymetry_individual - repmat(mu_columvol_ass_controls,[size(columvol_assymetry_individual,1),1]))./repmat(std_columvol_ass_controls,[size(columvol_assymetry_individual,1),1]);

% Assymetry T2 Signal
t2_signal_assymetry_individual        = 2*(t2_signal_individual(:,1:size(t2_signal_individual,2)/2) - t2_signal_individual(:,size(t2_signal_individual,2)/2+1:end))./(t2_signal_individual(:,1:size(t2_signal_individual,2)/2) + t2_signal_individual(:,size(t2_signal_individual,2)/2+1:end));
z_t2_signal_ass_patients_individual   = (t2_signal_assymetry_individual - repmat(mu_t2_signal_ass_controls,[size(t2_signal_assymetry_individual,1),1]))./repmat(std_t2_signal_ass_controls,[size(t2_signal_assymetry_individual,1),1]);

% Assymetry Ratio NUC Signal
ratio_nuc_signal_assymetry_individual         = 2*(ratio_nuc_signal_individual(:,1:size(ratio_nuc_signal_individual,2)/2) - ratio_nuc_signal_individual(:,size(ratio_nuc_signal_individual,2)/2+1:end))./(ratio_nuc_signal_individual(:,1:size(ratio_nuc_signal_individual,2)/2) + ratio_nuc_signal_individual(:,size(ratio_nuc_signal_individual,2)/2+1:end));
z_ratio_nuc_signal_ass_patients_individual    = (ratio_nuc_signal_assymetry_individual - repmat(mu_ratio_nuc_signal_ass_controls,[size(ratio_nuc_signal_assymetry_individual,1),1]))./repmat(std_ratio_nuc_signal_ass_controls,[size(ratio_nuc_signal_assymetry_individual,1),1]);

%% Run the Lateralization

% Build the input for the lateralizer
database_colvol_individual  = mean(z_columvol_ass_patients_individual(roi_colvol_final));
database_t2_individual      = mean(z_t2_signal_ass_patients_individual(roi_t2_final));
database_ratio_individual   = mean(z_ratio_nuc_signal_ass_patients_individual(roi_ratio_final));
database_t2_ratio(:,1)      = mean(z_t2_signal_ass_patients_individual(roi_t2_ratio_final));
database_t2_ratio(:,2)      = mean(z_ratio_nuc_signal_ass_patients_individual(roi_t2_ratio_final));

% Make predictions
[finalClassify_colvol_individual,tmp_posterior] = predict(lda_colvol_final,database_colvol_individual);
posterior_colvol_individual = tmp_posterior(2);

[finalClassify_t2_individual,tmp_posterior] = predict(lda_t2_final,database_t2_individual);
posterior_t2_individual = tmp_posterior(2);

[finalClassify_ratio_individual,tmp_posterior] = predict(lda_ratio_final,database_ratio_individual);
posterior_ratio_individual = tmp_posterior(2);

[finalClassify_t2_ratio_individual,tmp_posterior] = predict(lda_t2_ratio_final,database_t2_ratio_individual);
posterior_t2_ratio_individual = tmp_posterior(2);





