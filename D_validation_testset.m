%% Author: Benoit Caldairou, PhD
%% Mail: benoit.caldairou@mcgill.ca

%% Get Prisma Demographics
filename = '/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/06_FinalLateralization/Code_20210322_Currated/demographics_mni.txt';
delimiter = ',';

% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
Demographics_Test = [dataArray{1:end-1}];

% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Extract a few variables from demographics

prefix_test                 = Demographics_Test(:,1);
ids_test                    = Demographics_Test(:,2);
lateralization_truth_test   = Demographics_Test(:,3);
hs_clinical_test            = Demographics_Test(:,4);

%% Set feature directories
feature_dir='/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/ColumnVolume';
image_dir='';

%% Read the data (Columnar Volume)

ca_columvol_left_file_test   = strcat(featureDir,'/',prefix_test,'_',ids_test,'_L_CA_ColVol.txt');
ca_columvol_right_file_test  = strcat(featureDir,'/',prefix_test,'_',ids_test,'_R_CA_ColVol.txt');
dg_columvol_left_file_test   = strcat(featureDir,'/',prefix_test,'_',ids_test,'_L_DG_ColVol.txt');
dg_columvol_right_file_test  = strcat(featureDir,'/',prefix_test,'_',ids_test,'_R_DG_ColVol.txt');
sub_columvol_left_file_test  = strcat(featureDir,'/',prefix_test,'_',ids_test,'_L_SUB_ColVol.txt');
sub_columvol_right_file_test = strcat(featureDir,'/',prefix_test,'_',ids_test,'_R_SUB_ColVol.txt');

ca_columvol_left_test   = zeros(length(ids_test),10242);
ca_columvol_right_test  = zeros(length(ids_test),10242);
dg_columvol_left_test   = zeros(length(ids_test),5762);
dg_columvol_right_test  = zeros(length(ids_test),5762);
sub_columvol_left_test  = zeros(length(ids_test),5762);
sub_columvol_right_test = zeros(length(ids_test),5762);

for i=1:length(ids_test)

    ca_columvol_left_test(i,:) = SurfStatReadData(ca_columvol_left_file_test{i});
    ca_columvol_right_test(i,:) = SurfStatReadData(ca_columvol_right_file_test{i});
    dg_columvol_left_test(i,:) = SurfStatReadData(dg_columvol_left_file_test{i});
    dg_columvol_right_test(i,:) = SurfStatReadData(dg_columvol_right_file_test{i});
    sub_columvol_left_test(i,:) = SurfStatReadData(sub_columvol_left_file_test{i});
    sub_columvol_right_test(i,:) = SurfStatReadData(sub_columvol_right_file_test{i});

end

%% Read and normalize the data (T2 signal)

inDir_t2Signal  = '/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/T2Signal_Nuc';
inDir_t2Norm    = '/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/Ventriclemasks_Nuc';

ca_t2_left_file_test   = strcat(inDir_t2Signal,'/',prefix_test,'_',ids_test,'_L_CA_nnt2.txt');
ca_t2_right_file_test  = strcat(inDir_t2Signal,'/',prefix_test,'_',ids_test,'_R_CA_nnt2.txt');
dg_t2_left_file_test   = strcat(inDir_t2Signal,'/',prefix_test,'_',ids_test,'_L_DG_nnt2.txt');
dg_t2_right_file_test  = strcat(inDir_t2Signal,'/',prefix_test,'_',ids_test,'_R_DG_nnt2.txt');
sub_t2_left_file_test  = strcat(inDir_t2Signal,'/',prefix_test,'_',ids_test,'_L_SUB_nnt2.txt');
sub_t2_right_file_test = strcat(inDir_t2Signal,'/',prefix_test,'_',ids_test,'_R_SUB_nnt2.txt');
norm_t2_file_test      = strcat(inDir_t2Norm,'/',ids_test,'_T2norm.mat');

ca_t2_left_test   = zeros(length(ids_test),10242);
ca_t2_right_test  = zeros(length(ids_test),10242);
dg_t2_left_test   = zeros(length(ids_test),5762);
dg_t2_right_test  = zeros(length(ids_test),5762);
sub_t2_left_test  = zeros(length(ids_test),5762);
sub_t2_right_test = zeros(length(ids_test),5762);

for i=1:length(ids_test)

    ca_t2_left_test(i,:)   = SurfStatReadData(ca_t2_left_file_test{i});
    ca_t2_right_test(i,:)  = SurfStatReadData(ca_t2_right_file_test{i});
    dg_t2_left_test(i,:)   = SurfStatReadData(dg_t2_left_file_test{i});
    dg_t2_right_test(i,:)  = SurfStatReadData(dg_t2_right_file_test{i});
    sub_t2_left_test(i,:)  = SurfStatReadData(sub_t2_left_file_test{i});
    sub_t2_right_test(i,:) = SurfStatReadData(sub_t2_right_file_test{i});
    
    tmp = load(norm_t2_file_test{i});
    
    ca_t2_left_test(i,:)   = ca_t2_left_test(i,:)./tmp.Int_ref;
    ca_t2_right_test(i,:)  = ca_t2_right_test(i,:)./tmp.Int_ref;
    dg_t2_left_test(i,:)   = dg_t2_left_test(i,:)./tmp.Int_ref;
    dg_t2_right_test(i,:)  = dg_t2_right_test(i,:)./tmp.Int_ref;
    sub_t2_left_test(i,:)  = sub_t2_left_test(i,:)./tmp.Int_ref;
    sub_t2_right_test(i,:) = sub_t2_right_test(i,:)./tmp.Int_ref;

end

%% Read and normalize the ratio NUC data
inDir_ratioSignal = '/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/T2T1Ratio_Tal_Nuc';

ca_ratio_nuc_left_file_test   = strcat(inDir_ratioSignal,'/',prefix_test,'_',ids_test,'_L_CA_t2wt1wratio.txt');
ca_ratio_nuc_right_file_test  = strcat(inDir_ratioSignal,'/',prefix_test,'_',ids_test,'_R_CA_t2wt1wratio.txt');
dg_ratio_nuc_left_file_test   = strcat(inDir_ratioSignal,'/',prefix_test,'_',ids_test,'_L_DG_t2wt1wratio.txt');
dg_ratio_nuc_right_file_test  = strcat(inDir_ratioSignal,'/',prefix_test,'_',ids_test,'_R_DG_t2wt1wratio.txt');
sub_ratio_nuc_left_file_test  = strcat(inDir_ratioSignal,'/',prefix_test,'_',ids_test,'_L_SUB_t2wt1wratio.txt');
sub_ratio_nuc_right_file_test = strcat(inDir_ratioSignal,'/',prefix_test,'_',ids_test,'_R_SUB_t2wt1wratio.txt');

ca_ratio_nuc_left_test   = zeros(length(ids_test),10242);
ca_ratio_nuc_right_test  = zeros(length(ids_test),10242);
dg_ratio_nuc_left_test   = zeros(length(ids_test),5762);
dg_ratio_nuc_right_test  = zeros(length(ids_test),5762);
sub_ratio_nuc_left_test  = zeros(length(ids_test),5762);
sub_ratio_nuc_right_test = zeros(length(ids_test),5762);

for i=1:length(ids_test)

    ca_ratio_nuc_left_test(i,:)   = SurfStatReadData(ca_ratio_nuc_left_file_test{i});
    ca_ratio_nuc_right_test(i,:)  = SurfStatReadData(ca_ratio_nuc_right_file_test{i});
    dg_ratio_nuc_left_test(i,:)   = SurfStatReadData(dg_ratio_nuc_left_file_test{i});
    dg_ratio_nuc_right_test(i,:)  = SurfStatReadData(dg_ratio_nuc_right_file_test{i});
    sub_ratio_nuc_left_test(i,:)  = SurfStatReadData(sub_ratio_nuc_left_file_test{i});
    sub_ratio_nuc_right_test(i,:) = SurfStatReadData(sub_ratio_nuc_right_file_test{i});

end

%% z-score with respect to controls

% Global local volumes
columvol_test               = [ca_columvol_left_test,sub_columvol_left_test,dg_columvol_left_test,ca_columvol_right_test,sub_columvol_right_test,dg_columvol_right_test];
columvol_test               = SurfStatSmooth(columvol_test, template, 3);
z_columvol_patients_test    = (columvol_test - repmat(mu_columvol_controls,[size(columvol_test,1),1]))./repmat(std_columvol_controls,[size(columvol_test,1),1]);

% T2 signal
t2_signal_test              = [ca_t2_left_test,sub_t2_left_test,dg_t2_left_test,ca_t2_right_test,sub_t2_right_test,dg_t2_right_test];
t2_signal_test              = SurfStatSmooth(t2_signal_test, template, 3);
z_t2_signal_patients_test   = (t2_signal_test - repmat(mu_t2_signal_controls,[size(t2_signal_test,1),1]))./repmat(std_t2_signal_controls,[size(t2_signal_test,1),1]);

% Ratio NUC Signal
ratio_nuc_signal_test               = [ca_ratio_nuc_left_test,sub_ratio_nuc_left_test,dg_ratio_nuc_left_test,ca_ratio_nuc_right_test,sub_ratio_nuc_right_test,dg_ratio_nuc_right_test];
ratio_nuc_signal_test               = SurfStatSmooth(ratio_nuc_signal_test, template, 3);
z_ratio_nuc_signal_patients_test    = (ratio_nuc_signal_test - repmat(mu_ratio_nuc_signal_controls,[size(ratio_nuc_signal_test,1),1]))./repmat(std_ratio_nuc_signal_controls,[size(ratio_nuc_signal_test,1),1]);

% Assymetry columnar volumes
columvol_assymetry_test         = 2*(columvol_test(:,1:size(columvol_test,2)/2) - columvol_test(:,size(columvol_test,2)/2+1:end))./(columvol_test(:,1:size(columvol_test,2)/2) + columvol_test(:,size(columvol_test,2)/2+1:end));
z_columvol_ass_patients_test    = (columvol_assymetry_test - repmat(mu_columvol_ass_controls,[size(columvol_assymetry_test,1),1]))./repmat(std_columvol_ass_controls,[size(columvol_assymetry_test,1),1]);

% Assymetry T2 Signal
t2_signal_assymetry_test        = 2*(t2_signal_test(:,1:size(t2_signal_test,2)/2) - t2_signal_test(:,size(t2_signal_test,2)/2+1:end))./(t2_signal_test(:,1:size(t2_signal_test,2)/2) + t2_signal_test(:,size(t2_signal_test,2)/2+1:end));
z_t2_signal_ass_patients_test   = (t2_signal_assymetry_test - repmat(mu_t2_signal_ass_controls,[size(t2_signal_assymetry_test,1),1]))./repmat(std_t2_signal_ass_controls,[size(t2_signal_assymetry_test,1),1]);

% Assymetry Ratio NUC Signal
ratio_nuc_signal_assymetry_test         = 2*(ratio_nuc_signal_test(:,1:size(ratio_nuc_signal_test,2)/2) - ratio_nuc_signal_test(:,size(ratio_nuc_signal_test,2)/2+1:end))./(ratio_nuc_signal_test(:,1:size(ratio_nuc_signal_test,2)/2) + ratio_nuc_signal_test(:,size(ratio_nuc_signal_test,2)/2+1:end));
z_ratio_nuc_signal_ass_patients_test    = (ratio_nuc_signal_assymetry_test - repmat(mu_ratio_nuc_signal_ass_controls,[size(ratio_nuc_signal_assymetry_test,1),1]))./repmat(std_ratio_nuc_signal_ass_controls,[size(ratio_nuc_signal_assymetry_test,1),1]);

%% Switch from left right to ipsi contra

z_columvol_patients_ipsi_contra_test                = zeros(size(z_columvol_patients_test));
z_t2_signal_patients_ipsi_contra_test               = zeros(size(z_t2_signal_patients_test));
z_ratio_nuc_signal_patients_ipsi_contra_test        = zeros(size(z_ratio_nuc_signal_patients_test));

z_columvol_ass_patients_ipsi_contra_test            = zeros(size(z_columvol_ass_patients_test));
z_t2_signal_ass_patients_ipsi_contra_test           = zeros(size(z_t2_signal_ass_patients_test));
z_ratio_nuc_signal_ass_patients_ipsi_contra_test    = zeros(size(z_ratio_nuc_signal_ass_patients_test));

for i = 1:size(z_columvol_patients_test,1)

    if strcmp(lateralization_truth_test{i},'LTLE')
    
        % z score of raw values
        z_columvol_patients_ipsi_contra_test(i,:)           = z_columvol_patients_test(i,:);
        z_t2_signal_patients_ipsi_contra_test(i,:)          = z_t2_signal_patients_test(i,:);
        z_ratio_nuc_signal_patients_ipsi_contra_test(i,:)   = z_ratio_nuc_signal_patients_test(i,:);
        
        % z score of assymetries
        z_columvol_ass_patients_ipsi_contra_test(i,:)           = z_columvol_ass_patients_test(i,:);
        z_t2_signal_ass_patients_ipsi_contra_test(i,:)          = z_t2_signal_ass_patients_test(i,:);
        z_ratio_nuc_signal_ass_patients_ipsi_contra_test(i,:)   = z_ratio_nuc_signal_ass_patients_test(i,:);
        
    else
        
        % z score of raw values
        z_columvol_patients_ipsi_contra_test(i,1:size(z_columvol_patients_test,2)/2)        = z_columvol_patients_test(i,size(z_columvol_patients_test,2)/2+1:end);
        z_columvol_patients_ipsi_contra_test(i,size(z_columvol_patients_test,2)/2+1:end)    = z_columvol_patients_test(i,1:size(z_columvol_patients_test,2)/2);
        
        z_t2_signal_patients_ipsi_contra_test(i,1:size(z_t2_signal_patients_test,2)/2)      = z_t2_signal_patients_test(i,size(z_t2_signal_patients_test,2)/2+1:end);
        z_t2_signal_patients_ipsi_contra_test(i,size(z_t2_signal_patients_test,2)/2+1:end)  = z_t2_signal_patients_test(i,1:size(z_t2_signal_patients_test,2)/2);
        
        z_ratio_nuc_signal_patients_ipsi_contra_test(i,1:size(z_ratio_nuc_signal_patients_test,2)/2)        = z_ratio_nuc_signal_patients_test(i,size(z_ratio_nuc_signal_patients_test,2)/2+1:end);
        z_ratio_nuc_signal_patients_ipsi_contra_test(i,size(z_ratio_nuc_signal_patients_test,2)/2+1:end)    = z_ratio_nuc_signal_patients_test(i,1:size(z_ratio_nuc_signal_patients_test,2)/2);
        
        % z score of assymetries
        z_columvol_ass_patients_ipsi_contra_test(i,:)           = -z_columvol_ass_patients_test(i,:);
        z_t2_signal_ass_patients_ipsi_contra_test(i,:)          = -z_t2_signal_ass_patients_test(i,:);
        z_ratio_nuc_signal_ass_patients_ipsi_contra_test(i,:)   = -z_ratio_nuc_signal_ass_patients_test(i,:);
        
    end
    
end

%% Lateralization

testNumber = 3;
type = 'diagquadratic';

lateralization_truth_bin_test = strcmp(lateralization_truth_test,'LTLE');

roi_spam_colvol_test    = zeros(testNumber,21766);
roi_spam_t2_test        = zeros(testNumber,21766);
roi_spam_ratio_nuc_test = zeros(testNumber,21766);
roi_spam_t2_ratio_test  = zeros(testNumber,21766);

lda_colvol_test         = cell(testNumber,1);
lda_t2_test             = cell(testNumber,1);
lda_ratio_nuc_test      = cell(testNumber,1);
lda_t2_ratio_test       = cell(testNumber,1);

% Raw data labels
finalClassify_colvol_test      = zeros(length(lateralization_truth_bin_test),testNumber);
posterior_colvol_test          = zeros(length(lateralization_truth_bin_test),testNumber);

finalClassify_t2_test          = zeros(length(lateralization_truth_bin_test),testNumber);
posterior_t2_test              = zeros(length(lateralization_truth_bin_test),testNumber);

finalClassify_ratio_nuc_test   = zeros(length(lateralization_truth_bin_test),testNumber);
posterior_ratio_nuc_test       = zeros(length(lateralization_truth_bin_test),testNumber);

finalClassify_t2_ratio_test    = zeros(length(lateralization_truth_bin_test),testNumber);
posterior_t2_ratio_test        = zeros(length(lateralization_truth_bin_test),testNumber);


for i = 1:testNumber
    
    disp(i);
    partition = cvpartition(lateralization_bin_training,'HoldOut',0.2);
    
    % UNIMODAL SURFACE-BASED
    % ------------------------- ColVol ------------------------------ %
    % t-test to get t-maps beween ipsi and contra for training set
    [~,~,~,stats_nohs] = ttest2(z_columvol_patients_ipsi_contra(partition.training,1:size(z_columvol_patients,2)/2),z_columvol_patients_ipsi_contra(partition.training,size(z_columvol_patients,2)/2+1:end), 0.025, 'both','unequal');
    tmap = abs(stats_nohs.tstat);
    
    % Looks for optimal t in training set
    training_data  = z_columvol_ass_patients(partition.training,:);
    training_group = lateralization_bin_training(partition.training,:);
    t_optimal = optimal_thres_lookup(tmap,training_data,training_group);
    
    % Build data for this specific t
    roi = tmap >= t_optimal;
    roi_spam_colvol_test(i,:) = roi;
    database_train = mean(z_columvol_ass_patients(partition.training,roi),2);
    database_test = mean(z_columvol_ass_patients_test(:,roi),2);
    
    % Build the model and classify
    lda_colvol_test{i} = fitcdiscr(database_train,lateralization_bin_training(partition.training),'DiscrimType',type);
    [finalClassify_colvol_test(:,i),tmp_posterior] = predict(lda_colvol_test{i},database_test);
    posterior_colvol_test(:,i) = tmp_posterior(:,2);
    
    % --------------------------- T2 ------------------------------- %
    % t-test to get t-maps beween ipsi and contra for training set
    [~,~,~,stats_nohs] = ttest2(z_t2_signal_patients_ipsi_contra(partition.training,1:size(z_t2_signal_patients,2)/2),z_t2_signal_patients_ipsi_contra(partition.training,size(z_t2_signal_patients,2)/2+1:end), 0.025, 'both','unequal');
    tmap = abs(stats_nohs.tstat);
    
    % Looks for optimal t in training set
    training_data  = z_t2_signal_ass_patients(partition.training,:);
    training_group = lateralization_bin_training(partition.training);
    t_optimal = optimal_thres_lookup(tmap,training_data,training_group);
    
    % Build data for this specific t
    roi = tmap >= t_optimal;
    database_train = mean(z_t2_signal_ass_patients(partition.training,roi),2);
    database_test = mean(z_t2_signal_ass_patients_test(:,roi),2);
    
    % Build the model and classify
    lda_t2_test{i} = fitcdiscr(database_train,lateralization_bin_training(partition.training),'DiscrimType',type);
    [finalClassify_t2_test(:,i),tmp_posterior] = predict(lda_t2_test{i},database_test);
    posterior_t2_test(:,i) = tmp_posterior(:,2);
    
    % --------------------------- RATIO NUC ----------------------------- %
    % t-test to get t-maps beween ipsi and contra for training set
    [~,~,~,stats_nohs] = ttest2(z_ratio_nuc_signal_patients(partition.training,1:size(z_ratio_nuc_signal_patients,2)/2),z_ratio_nuc_signal_patients_ipsi_contra(partition.training,size(z_ratio_nuc_signal_patients,2)/2+1:end), 0.025, 'both','unequal');
    tmap = abs(stats_nohs.tstat);
    
    % Looks for optimal t in training set
    training_data  = z_ratio_nuc_signal_ass_patients(partition.training,:);
    training_group = lateralization_bin_training(partition.training);
    t_optimal = optimal_thres_lookup(tmap,training_data,training_group);
    
    % Build data for this specific t
    roi = tmap >= t_optimal;
    database_train = mean(z_ratio_nuc_signal_ass_patients(partition.training,roi),2);
    database_test = mean(z_ratio_nuc_signal_ass_patients_test(:,roi),2);
    
    % Build the model and classify
    lda_ratio_nuc_test{i} = fitcdiscr(database_train,lateralization_bin_training(partition.training),'DiscrimType',type);
    [finalClassify_ratio_nuc_test(:,i),tmp_posterior] = predict(lda_ratio_nuc_test{i},database_test);
    posterior_ratio_nuc_test(:,i) = tmp_posterior(:,2);
        
    % ------------- MULTIVARIATE T2 RATIO -------------------------- %
    clear forT2test_z_ipsi_contra_columnarvolume_training forT2test_z_ipsi_contra_ratiosignal_training forT2test_z_ipsi_contra_all_training;
    
    % t-test to get t-maps beween ipsi and contra for training set
    forT2test_z_ipsi_contra_t2signal_training = [z_t2_signal_patients_ipsi_contra(partition.training,1:size(z_t2_signal_patients,2)/2);z_t2_signal_patients_ipsi_contra(partition.training,size(z_t2_signal_patients,2)/2+1:end)];
    forT2test_z_ipsi_contra_ratiosignal_training = [z_ratio_nuc_signal_patients_ipsi_contra(partition.training,1:size(z_ratio_nuc_signal_patients,2)/2);z_ratio_nuc_signal_patients_ipsi_contra(partition.training,size(z_ratio_nuc_signal_patients,2)/2+1:end)];
    forT2test_z_ipsi_contra_all_training = cat(3,forT2test_z_ipsi_contra_t2signal_training,forT2test_z_ipsi_contra_ratiosignal_training);
    
    groups_ipsi_contra = {};
    groups_ipsi_contra(1:size(forT2test_z_ipsi_contra_t2signal_training,1)/2) = {'IPSI'};
    groups_ipsi_contra(size(forT2test_z_ipsi_contra_t2signal_training,1)/2+1:size(forT2test_z_ipsi_contra_t2signal_training,1)) = {'CONTRA'};
    GROUPS_IPSI_CONTRA = term(cellstr(groups_ipsi_contra));
    
    % Put a very simple model
    M = 1 + GROUPS_IPSI_CONTRA;
    % Hotelling T2 to get T-maps beween ipsi and contra
    slm_all = SurfStatLinMod(forT2test_z_ipsi_contra_all_training, M, template_uni);
    slm_all = SurfStatT(slm_all,GROUPS_IPSI_CONTRA.IPSI-GROUPS_IPSI_CONTRA.CONTRA);% A few more variables
    
    tmap = abs(slm_all.t);
    
    clear training_data database database_prisma;
    
    % Looks for optimal t in training set
    training_data(:,:,1)  = z_t2_signal_ass_patients(partition.training,:);
    training_data(:,:,2)  = z_ratio_nuc_signal_ass_patients(partition.training,:);
    training_group = lateralization_bin_training(partition.training);
    t_optimal = optimal_thres_lookup(tmap,training_data,training_group);
    
    % Build data for this specific t
    roi = tmap >= t_optimal;
    database_train(:,1) = mean(z_t2_signal_ass_patients(partition.training,roi),2);
    database_train(:,2) = mean(z_ratio_nuc_signal_ass_patients(partition.training,roi),2);
    database_test(:,1) = mean(z_t2_signal_ass_patients_test(:,roi),2);
    database_test(:,2) = mean(z_ratio_nuc_signal_ass_patients_test(:,roi),2);
    
    % Build the model and classify
    lda_t2_ratio_test{i} = fitcdiscr(database_train,lateralization_bin_training(partition.training),'DiscrimType','diagquadratic');
    [finalClassify_t2_ratio_test(:,i),tmp_posterior] = predict(lda_t2_ratio_test{i},database_test);
    posterior_t2_ratio_test(:,i) = tmp_posterior(:,2);
    
end


%% Check the laterlization accuracy

nb_test = length(lateralization_truth_bin_test);

for i = 1:size(finalClassify_colvol_test,2)
    colvol_test_comparison(:,i) = finalClassify_colvol_test(:,i) == lateralization_truth_bin_test;
    t2_test_comparison(:,i) = finalClassify_t2_test(:,i) == lateralization_truth_bin_test;
    ratio_test_comparison(:,i) = finalClassify_ratio_nuc_test(:,i) == lateralization_truth_bin_test;
    t2_ratio_test_comparison(:,i) = finalClassify_t2_ratio_test(:,i) == lateralization_truth_bin_test;
end

disp('Individual Performance')
["IDs", "Truth", "ColVol", "T2", "Ratio Nuc", "T2 Ratio";
ids_test, lateralization_truth_test, num2cell(sum(colvol_test_comparison,2)), num2cell(sum(t2_test_comparison,2)), num2cell(sum(ratio_test_comparison,2)), num2cell(sum(t2_ratio_test_comparison,2))]

[mean(sum(colvol_test_comparison)),std(sum(colvol_test_comparison)),mean(sum(colvol_test_comparison))/nb_test,std(sum(colvol_test_comparison))/nb_test]
[mean(sum(t2_test_comparison)),std(sum(t2_test_comparison)),mean(sum(t2_test_comparison))/nb_test,std(sum(t2_test_comparison))/nb_test]
[mean(sum(ratio_test_comparison)),std(sum(ratio_test_comparison)),mean(sum(ratio_test_comparison))/nb_test,std(sum(ratio_test_comparison))/nb_test]
[mean(sum(t2_ratio_test_comparison)),std(sum(t2_ratio_test_comparison)),mean(sum(t2_ratio_test_comparison))/nb_test,std(sum(t2_ratio_test_comparison))/nb_test]

nb_test = length(lateralization_truth_bin_test(strcmp(hs_clinical_test,"yes")));
[mean(sum(colvol_test_comparison(strcmp(hs_clinical_test,"yes"),:))),std(sum(colvol_test_comparison(strcmp(hs_clinical_test,"yes"),:))),mean(sum(colvol_test_comparison(strcmp(hs_clinical_test,"yes"),:)))/nb_test,std(sum(colvol_test_comparison(strcmp(hs_clinical_test,"yes"),:)))/nb_test]
[mean(sum(t2_test_comparison(strcmp(hs_clinical_test,"yes"),:))),std(sum(t2_test_comparison(strcmp(hs_clinical_test,"yes"),:))),mean(sum(t2_test_comparison(strcmp(hs_clinical_test,"yes"),:)))/nb_test,std(sum(t2_test_comparison(strcmp(hs_clinical_test,"yes"),:)))/nb_test]
[mean(sum(ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:))),std(sum(ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:))),mean(sum(ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:)))/nb_test,std(sum(ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:)))/nb_test]
[mean(sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:))),std(sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:))),mean(sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:)))/nb_test,std(sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:)))/nb_test]

nb_test = length(lateralization_truth_bin_test(strcmp(hs_clinical_test,"no")));
[mean(sum(colvol_test_comparison(strcmp(hs_clinical_test,"no"),:))),std(sum(colvol_test_comparison(strcmp(hs_clinical_test,"no"),:))),mean(sum(colvol_test_comparison(strcmp(hs_clinical_test,"no"),:)))/nb_test,std(sum(colvol_test_comparison(strcmp(hs_clinical_test,"no"),:)))/nb_test]
[mean(sum(t2_test_comparison(strcmp(hs_clinical_test,"no"),:))),std(sum(t2_test_comparison(strcmp(hs_clinical_test,"no"),:))),mean(sum(t2_test_comparison(strcmp(hs_clinical_test,"no"),:)))/nb_test,std(sum(t2_test_comparison(strcmp(hs_clinical_test,"no"),:)))/nb_test]
[mean(sum(ratio_test_comparison(strcmp(hs_clinical_test,"no"),:))),std(sum(ratio_test_comparison(strcmp(hs_clinical_test,"no"),:))),mean(sum(ratio_test_comparison(strcmp(hs_clinical_test,"no"),:)))/nb_test,std(sum(ratio_test_comparison(strcmp(hs_clinical_test,"no"),:)))/nb_test]
[mean(sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"no"),:))),std(sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"no"),:))),mean(sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"no"),:)))/nb_test,std(sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"no"),:)))/nb_test]

%% Significance Test
disp('------- All Patients -------');
disp('Friedman Test');
all_results = [ sum(colvol_test_comparison)',... 
                sum(t2_test_comparison)',... 
                sum(ratio_test_comparison)',...
                sum(t2_ratio_test_comparison)',...
                ];
            
[p_all_friedman,tbl,stats_all_friedman] = friedman(all_results);
c_all_friedman = multcompare(stats_all_friedman,'CType','bonferroni');

disp('--------- MRI+ Patients ---------');
all_results = [ sum(colvol_test_comparison(strcmp(hs_clinical_test,"yes"),:))',... 
                sum(t2_test_comparison(strcmp(hs_clinical_test,"yes"),:))',... 
                sum(ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:))',...
                sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"yes"),:))',...
                ];
            
[p_all_friedman,tbl,stats_all_friedman] = friedman(all_results);
c_all_friedman = multcompare(stats_all_friedman,'CType','bonferroni');

disp('--------- MRI- Patients ---------');
all_results = [ sum(colvol_test_comparison(strcmp(hs_clinical_test,"no"),:))',... 
                sum(t2_test_comparison(strcmp(hs_clinical_test,"no"),:))',... 
                sum(ratio_test_comparison(strcmp(hs_clinical_test,"no"),:))',...
                sum(t2_ratio_test_comparison(strcmp(hs_clinical_test,"no"),:))',...
                ];
            
[p_all_friedman,tbl,stats_all_friedman] = friedman(all_results);
c_all_friedman = multcompare(stats_all_friedman,'CType','bonferroni');


%% ROC curves
figure(2);

meanLineWidth = 4.0;
meanLineColor = '#ff0000';
meanLineStyle = '-';

individualLineColor = '#0080ff';
individualLineStyle = ':';

textFontSize = 24;

nb_ltle = sum(lateralization_truth_bin_test);
nb_rtle = sum(~lateralization_truth_bin_test);

nb_ltle_hs      = sum(lateralization_truth_bin_test & strcmp(hs_clinical_test,'yes'));
nb_ltle_nonhs   = sum(lateralization_truth_bin_test & strcmp(hs_clinical_test,'no'));

nb_rtle_hs      = sum(~lateralization_truth_bin_test & strcmp(hs_clinical_test,'yes'));
nb_rtle_nonhs   = sum(~lateralization_truth_bin_test & strcmp(hs_clinical_test,'no'));

% Colvol Subplot

clear -regex ltle_positive_* rtle_positive_* auc* ;

ax= subplot(2,2,1,'FontSize',textFontSize); hold on;
ax.XLabel.String = 'LTLE FPR'; ax.YLabel.String = 'LTLE TPR';

plot(0:0.001:1,0:0.001:1,'-k');

for repeat = 1:testNumber

    tmp_posterior = posterior_colvol_test(:,repeat);
    t=0;
    
    for thres=0:0.001:1
    
        t = t+1;
        tmp_classify = tmp_posterior > thres;
        tmp_result = lateralization_truth_bin_test == tmp_classify;
        
        tmp_result_hs = tmp_result(strcmp(hs_clinical_test,'yes'));
        tmp_result_nonhs = tmp_result(strcmp(hs_clinical_test,'no'));
        
        ltle_positive_colvol_test(repeat,t) = sum(tmp_result & lateralization_truth_bin_test);
        rtle_positive_colvol_test(repeat,t) = sum(tmp_result & ~lateralization_truth_bin_test);
        
        ltle_positive_colvol_test_hs(repeat,t) = sum(tmp_result_hs & lateralization_truth_bin_test(strcmp(hs_clinical_test,'yes')));
        rtle_positive_colvol_test_hs(repeat,t) = sum(tmp_result_hs & ~lateralization_truth_bin_test(strcmp(hs_clinical_test,'yes')));
        
        ltle_positive_colvol_test_nonhs(repeat,t) = sum(tmp_result_nonhs & lateralization_truth_bin_test(strcmp(hs_clinical_test,'no')));
        rtle_positive_colvol_test_nonhs(repeat,t) = sum(tmp_result_nonhs & ~lateralization_truth_bin_test(strcmp(hs_clinical_test,'no')));
        
        %ltle_positive_colvol_operated_hs(repeat,t) = ltle
    end

    ltle_positive_colvol_test(repeat,:) = ltle_positive_colvol_test(repeat,:)/nb_ltle;
    rtle_positive_colvol_test(repeat,:) = rtle_positive_colvol_test(repeat,:)/nb_rtle;
    
    ltle_positive_colvol_test_hs(repeat,:) = ltle_positive_colvol_test_hs(repeat,:)/nb_ltle_hs;
    rtle_positive_colvol_test_hs(repeat,:) = rtle_positive_colvol_test_hs(repeat,:)/nb_rtle_hs;
    
    ltle_positive_colvol_test_nonhs(repeat,:) = ltle_positive_colvol_test_nonhs(repeat,:)/nb_ltle_nonhs;
    rtle_positive_colvol_test_nonhs(repeat,:) = rtle_positive_colvol_test_nonhs(repeat,:)/nb_rtle_nonhs;
    
    auc_colvol_test(repeat)         = -trapz(1-rtle_positive_colvol_test(repeat,:),ltle_positive_colvol_test(repeat,:));
    auc_colvol_test_hs(repeat,:)    = -trapz(1-rtle_positive_colvol_test_hs(repeat,:),ltle_positive_colvol_test_hs(repeat,:));
    auc_colvol_test_nonhs(repeat,:) = -trapz(1-rtle_positive_colvol_test_nonhs(repeat,:),ltle_positive_colvol_test_nonhs(repeat,:));
    
    plot(1-rtle_positive_colvol_test_nonhs(repeat,:),ltle_positive_colvol_test_nonhs(repeat,:),'LineStyle',individualLineStyle,'Color', individualLineColor);

end

plot(mean(1-rtle_positive_colvol_test_nonhs),mean(ltle_positive_colvol_test_nonhs),'LineStyle',meanLineStyle,'Color',meanLineColor,'LineWidth',meanLineWidth);
text(0.35,0.1,['AUC = ', num2str(mean(auc_colvol_test_nonhs),'%0.2g')],'FontSize', textFontSize);

hold off;

% T2 Subplot

subplot(2,2,2,'FontSize',textFontSize); hold on;
plot(0:0.001:1,0:0.001:1,'-k');

for repeat = 1:testNumber

    tmp_posterior = posterior_t2_test(:,repeat);
    t=0;
    
    for thres=0:0.001:1
    
        t = t+1;
        tmp_classify = tmp_posterior > thres;
        tmp_result = lateralization_truth_bin_test == tmp_classify;
        
        tmp_result_hs = tmp_result(strcmp(hs_clinical_test,'yes'));
        tmp_result_nonhs = tmp_result(strcmp(hs_clinical_test,'no'));
        
        ltle_positive_t2_test(repeat,t) = sum(tmp_result & lateralization_truth_bin_test);
        rtle_positive_t2_test(repeat,t) = sum(tmp_result & ~lateralization_truth_bin_test);
        
        ltle_positive_t2_test_hs(repeat,t) = sum(tmp_result_hs & lateralization_truth_bin_test(strcmp(hs_clinical_test,'yes')));
        rtle_positive_t2_test_hs(repeat,t) = sum(tmp_result_hs & ~lateralization_truth_bin_test(strcmp(hs_clinical_test,'yes')));
        
        ltle_positive_t2_test_nonhs(repeat,t) = sum(tmp_result_nonhs & lateralization_truth_bin_test(strcmp(hs_clinical_test,'no')));
        rtle_positive_t2_test_nonhs(repeat,t) = sum(tmp_result_nonhs & ~lateralization_truth_bin_test(strcmp(hs_clinical_test,'no')));
        
    end

    ltle_positive_t2_test(repeat,:) = ltle_positive_t2_test(repeat,:)/nb_ltle;
    rtle_positive_t2_test(repeat,:) = rtle_positive_t2_test(repeat,:)/nb_rtle;
    
    ltle_positive_t2_test_hs(repeat,:) = ltle_positive_t2_test_hs(repeat,:)/nb_ltle_hs;
    rtle_positive_t2_test_hs(repeat,:) = rtle_positive_t2_test_hs(repeat,:)/nb_rtle_hs;
    
    ltle_positive_t2_test_nonhs(repeat,:) = ltle_positive_t2_test_nonhs(repeat,:)/nb_ltle_nonhs;
    rtle_positive_t2_test_nonhs(repeat,:) = rtle_positive_t2_test_nonhs(repeat,:)/nb_rtle_nonhs;
    
    auc_t2_test(repeat)         = -trapz(1-rtle_positive_t2_test(repeat,:),ltle_positive_t2_test(repeat,:));
    auc_t2_test_hs(repeat,:)    = -trapz(1-rtle_positive_t2_test_hs(repeat,:),ltle_positive_t2_test_hs(repeat,:));
    auc_t2_test_nonhs(repeat,:) = -trapz(1-rtle_positive_t2_test_nonhs(repeat,:),ltle_positive_t2_test_nonhs(repeat,:));
    
    plot(1-rtle_positive_t2_test_nonhs(repeat,:),ltle_positive_t2_test_nonhs(repeat,:),'LineStyle',individualLineStyle,'Color', individualLineColor);

end

plot(mean(1-rtle_positive_t2_test_nonhs),mean(ltle_positive_t2_test_nonhs),'LineStyle',meanLineStyle,'Color',meanLineColor,'LineWidth',meanLineWidth);
text(0.35,0.1,['AUC = ', num2str(mean(auc_t2_test_nonhs),'%0.2g')],'FontSize', textFontSize);

hold off;

% Ratio Subplot

subplot(2,2,3,'FontSize',textFontSize); hold on;
plot(0:0.001:1,0:0.001:1,'-k');

for repeat = 1:testNumber

    tmp_posterior = posterior_ratio_nuc_test(:,repeat);
    t=0;
    
    for thres=0:0.001:1
    
        t = t+1;
        tmp_classify = tmp_posterior > thres;
        tmp_result = lateralization_truth_bin_test == tmp_classify;
        
        tmp_result_hs = tmp_result(strcmp(hs_clinical_test,'yes'));
        tmp_result_nonhs = tmp_result(strcmp(hs_clinical_test,'no'));
        
        ltle_positive_ratio_nuc_test(repeat,t) = sum(tmp_result & lateralization_truth_bin_test);
        rtle_positive_ratio_nuc_test(repeat,t) = sum(tmp_result & ~lateralization_truth_bin_test);
        
        ltle_positive_ratio_nuc_test_hs(repeat,t) = sum(tmp_result_hs & lateralization_truth_bin_test(strcmp(hs_clinical_test,'yes')));
        rtle_positive_ratio_nuc_test_hs(repeat,t) = sum(tmp_result_hs & ~lateralization_truth_bin_test(strcmp(hs_clinical_test,'yes')));
        
        ltle_positive_ratio_nuc_test_nonhs(repeat,t) = sum(tmp_result_nonhs & lateralization_truth_bin_test(strcmp(hs_clinical_test,'no')));
        rtle_positive_ratio_nuc_test_nonhs(repeat,t) = sum(tmp_result_nonhs & ~lateralization_truth_bin_test(strcmp(hs_clinical_test,'no')));
        
        
    end

    ltle_positive_ratio_nuc_test(repeat,:) = ltle_positive_ratio_nuc_test(repeat,:)/nb_ltle;
    rtle_positive_ratio_nuc_test(repeat,:) = rtle_positive_ratio_nuc_test(repeat,:)/nb_rtle;
    
    ltle_positive_ratio_nuc_test_hs(repeat,:) = ltle_positive_ratio_nuc_test_hs(repeat,:)/nb_ltle_hs;
    rtle_positive_ratio_nuc_test_hs(repeat,:) = rtle_positive_ratio_nuc_test_hs(repeat,:)/nb_rtle_hs;
    
    ltle_positive_ratio_nuc_test_nonhs(repeat,:) = ltle_positive_ratio_nuc_test_nonhs(repeat,:)/nb_ltle_nonhs;
    rtle_positive_ratio_nuc_test_nonhs(repeat,:) = rtle_positive_ratio_nuc_test_nonhs(repeat,:)/nb_rtle_nonhs;
    
    auc_ratio_nuc_test(repeat)          = -trapz(1-rtle_positive_ratio_nuc_test(repeat,:),ltle_positive_ratio_nuc_test(repeat,:));
    auc_ratio_nuc_test_hs(repeat,:)     = -trapz(1-rtle_positive_ratio_nuc_test_hs(repeat,:),ltle_positive_ratio_nuc_test_hs(repeat,:));
    auc_ratio_nuc_test_nonhs(repeat,:)  = -trapz(1-rtle_positive_ratio_nuc_test_nonhs(repeat,:),ltle_positive_ratio_nuc_test_nonhs(repeat,:));
    
    plot(1-rtle_positive_ratio_nuc_test_nonhs(repeat,:),ltle_positive_ratio_nuc_test_nonhs(repeat,:),'LineStyle',individualLineStyle,'Color', individualLineColor);

end

plot(mean(1-rtle_positive_ratio_nuc_test_nonhs),mean(ltle_positive_ratio_nuc_test_nonhs),'LineStyle',meanLineStyle,'Color',meanLineColor,'LineWidth',meanLineWidth);
text(0.35,0.1,['AUC = ', num2str(mean(auc_ratio_nuc_test_nonhs),'%0.2g')],'FontSize', textFontSize);

hold off;

% T2/Ratio Subplot

subplot(2,2,4,'FontSize',textFontSize); hold on;
plot(0:0.001:1,0:0.001:1,'-k');

for repeat = 1:testNumber

    tmp_posterior = posterior_t2_ratio_test(:,repeat);
    t=0;
    
    for thres=0:0.001:1
    
        t = t+1;
        tmp_classify = tmp_posterior > thres;
        tmp_result = lateralization_truth_bin_test == tmp_classify;
        
        tmp_result_hs = tmp_result(strcmp(hs_clinical_test,'yes'));
        tmp_result_nonhs = tmp_result(strcmp(hs_clinical_test,'no'));
        
        ltle_positive_t2_ratio_test(repeat,t) = sum(tmp_result & lateralization_truth_bin_test);
        rtle_positive_t2_ratio_test(repeat,t) = sum(tmp_result & ~lateralization_truth_bin_test);
        
        ltle_positive_t2_ratio_test_hs(repeat,t) = sum(tmp_result_hs & lateralization_truth_bin_test(strcmp(hs_clinical_test,'yes')));
        rtle_positive_t2_ratio_test_hs(repeat,t) = sum(tmp_result_hs & ~lateralization_truth_bin_test(strcmp(hs_clinical_test,'yes')));
        
        ltle_positive_t2_ratio_test_nonhs(repeat,t) = sum(tmp_result_nonhs & lateralization_truth_bin_test(strcmp(hs_clinical_test,'no')));
        rtle_positive_t2_ratio_test_nonhs(repeat,t) = sum(tmp_result_nonhs & ~lateralization_truth_bin_test(strcmp(hs_clinical_test,'no')));
        
    end

    ltle_positive_t2_ratio_test(repeat,:) = ltle_positive_t2_ratio_test(repeat,:)/nb_ltle;
    rtle_positive_t2_ratio_test(repeat,:) = rtle_positive_t2_ratio_test(repeat,:)/nb_rtle;
    
    ltle_positive_t2_ratio_test_hs(repeat,:) = ltle_positive_t2_ratio_test_hs(repeat,:)/nb_ltle_hs;
    rtle_positive_t2_ratio_test_hs(repeat,:) = rtle_positive_t2_ratio_test_hs(repeat,:)/nb_rtle_hs;
    
    ltle_positive_t2_ratio_test_nonhs(repeat,:) = ltle_positive_t2_ratio_test_nonhs(repeat,:)/nb_ltle_nonhs;
    rtle_positive_t2_ratio_test_nonhs(repeat,:) = rtle_positive_t2_ratio_test_nonhs(repeat,:)/nb_rtle_nonhs;
    
    auc_t2_ratio_test(repeat)           = -trapz(1-rtle_positive_t2_ratio_test(repeat,:),ltle_positive_t2_ratio_test(repeat,:));
    auc_t2_ratio_test_hs(repeat,:)      = -trapz(1-rtle_positive_t2_ratio_test_hs(repeat,:),ltle_positive_t2_ratio_test_hs(repeat,:));
    auc_t2_ratio_test_nonhs(repeat,:)   = -trapz(1-rtle_positive_t2_ratio_test_nonhs(repeat,:),ltle_positive_t2_ratio_test_nonhs(repeat,:));
    
    plot(1-rtle_positive_t2_ratio_test_nonhs(repeat,:),ltle_positive_t2_ratio_test_nonhs(repeat,:),'LineStyle',individualLineStyle,'Color', individualLineColor);

end

plot(mean(1-rtle_positive_t2_ratio_test_nonhs),mean(ltle_positive_t2_ratio_test_nonhs),'LineStyle',meanLineStyle,'Color',meanLineColor,'LineWidth',meanLineWidth);
text(0.35,0.1,['AUC = ', num2str(mean(auc_t2_ratio_test_nonhs),'%0.2g')],'FontSize', textFontSize);

hold off;

