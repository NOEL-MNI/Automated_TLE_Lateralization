%% Author: Benoit Caldairou, PhD
%% Mail: benoit.caldairou@mcgill.ca

%% Import Demographics.
% Script for importing data from the following text file:
%
%    /host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/06_FinalLateralization/demographics_20190709.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2019/07/09 14:38:33

% Initialize variables.
filename = '/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/06_FinalLateralization/Code_20210322_Currated/demographics_training_set.txt';
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
Demographics = [dataArray{1:end-1}];

% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Set lib path
addpath('/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/06_FinalLateralization/Code_20210322_Currated/lib');

%% Extract a few variables
% prefix: prefix establishing population belonging
% ids: individual ids
% groups_training: LTLE or RTLE or Control
% hs_clinical_training: MRI+ ("yes") or MRI- ("no") or Control ("na")

prefix                  = Demographics(:,1);
ids                     = Demographics(:,2);
groups_training         = Demographics(:,3);
hs_clinical_training    = Demographics(:,4);


% Separate controls from patients
controls_group                  = strcmp('Controls',groups_training);
tle_group_training              = ~controls_group;
lateralization_truth_training   = groups_training(tle_group_training);
hs_patients_training            = hs_clinical_training(tle_group_training);
lateralization_bin_training     = logical(strcmp(lateralization_truth_training,'LTLE')); % Binarize patient group: 1 for LTLE, 0 for RTLE

ids_patients = ids(tle_group_training);

% Load surface template
template = SurfStatReadSurf1('template.obj');
template_uni.coord = template.coord(:,1:21766);
template_uni.normal = template.normal(:,1:21766);
template_uni.tri = template.tri(1:43520,:);
template_uni.colr = template.colr;

cohen_colormap = flipud(brewermap(1000,'RdBu'));

%% Image and Data directories
feature_dir='/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling';
image_dir='';

%% Read the Columnar Volume

ca_columvol_left_file   = strcat(feature_dir,'/',prefix,'_',ids,'_L_CA_ColVol.txt');
ca_columvol_right_file  = strcat(feature_dir,'/',prefix,'_',ids,'_R_CA_ColVol.txt');
dg_columvol_left_file   = strcat(feature_dir,'/',prefix,'_',ids,'_L_DG_ColVol.txt');
dg_columvol_right_file  = strcat(feature_dir,'/',prefix,'_',ids,'_R_DG_ColVol.txt');
sub_columvol_left_file  = strcat(feature_dir,'/',prefix,'_',ids,'_L_SUB_ColVol.txt');
sub_columvol_right_file = strcat(feature_dir,'/',prefix,'_',ids,'_R_SUB_ColVol.txt');

% Depends on the number of vertices. Might have to be adjusted.
ca_columvol_left   = zeros(length(ids),10242);
ca_columvol_right  = zeros(length(ids),10242);
dg_columvol_left   = zeros(length(ids),5762);
dg_columvol_right  = zeros(length(ids),5762);
sub_columvol_left  = zeros(length(ids),5762);
sub_columvol_right = zeros(length(ids),5762);

for i=1:length(ids)

    ca_columvol_left(i,:) = SurfStatReadData(ca_columvol_left_file{i});
    ca_columvol_right(i,:) = SurfStatReadData(ca_columvol_right_file{i});
    dg_columvol_left(i,:) = SurfStatReadData(dg_columvol_left_file{i});
    dg_columvol_right(i,:) = SurfStatReadData(dg_columvol_right_file{i});
    sub_columvol_left(i,:) = SurfStatReadData(sub_columvol_left_file{i});
    sub_columvol_right(i,:) = SurfStatReadData(sub_columvol_right_file{i});

end

%% Read and normalize the data (T2 signal)

ca_t2_left_file   = strcat(feature_dir,'/',prefix,'_',ids,'_L_CA_nnt2.txt');
ca_t2_right_file  = strcat(feature_dir,'/',prefix,'_',ids,'_R_CA_nnt2.txt');
dg_t2_left_file   = strcat(feature_dir,'/',prefix,'_',ids,'_L_DG_nnt2.txt');
dg_t2_right_file  = strcat(feature_dir,'/',prefix,'_',ids,'_R_DG_nnt2.txt');
sub_t2_left_file  = strcat(feature_dir,'/',prefix,'_',ids,'_L_SUB_nnt2.txt');
sub_t2_right_file = strcat(feature_dir,'/',prefix,'_',ids,'_R_SUB_nnt2.txt');
norm_t2_file      = strcat(image_dir,'/',ids,'_T2norm.mat');

ca_t2_left   = zeros(length(ids),10242);
ca_t2_right  = zeros(length(ids),10242);
dg_t2_left   = zeros(length(ids),5762);
dg_t2_right  = zeros(length(ids),5762);
sub_t2_left  = zeros(length(ids),5762);
sub_t2_right = zeros(length(ids),5762);

for i=1:length(ids)

    ca_t2_left(i,:)   = SurfStatReadData(ca_t2_left_file{i});
    ca_t2_right(i,:)  = SurfStatReadData(ca_t2_right_file{i});
    dg_t2_left(i,:)   = SurfStatReadData(dg_t2_left_file{i});
    dg_t2_right(i,:)  = SurfStatReadData(dg_t2_right_file{i});
    sub_t2_left(i,:)  = SurfStatReadData(sub_t2_left_file{i});
    sub_t2_right(i,:) = SurfStatReadData(sub_t2_right_file{i});
    
    tmp = load(norm_t2_file{i});
    
    ca_t2_left(i,:)   = ca_t2_left(i,:)./tmp.Int_ref;
    ca_t2_right(i,:)  = ca_t2_right(i,:)./tmp.Int_ref;
    dg_t2_left(i,:)   = dg_t2_left(i,:)./tmp.Int_ref;
    dg_t2_right(i,:)  = dg_t2_right(i,:)./tmp.Int_ref;
    sub_t2_left(i,:)  = sub_t2_left(i,:)./tmp.Int_ref;
    sub_t2_right(i,:) = sub_t2_right(i,:)./tmp.Int_ref;

end

%% Read flair/t1 Ratio after non-uniformity correction

ca_ratio_nuc_left_file   = strcat(feature_dir,'/',prefix,'_',ids,'_L_CA_t2wt1wratio.txt');
ca_ratio_nuc_right_file  = strcat(feature_dir,'/',prefix,'_',ids,'_R_CA_t2wt1wratio.txt');
dg_ratio_nuc_left_file   = strcat(feature_dir,'/',prefix,'_',ids,'_L_DG_t2wt1wratio.txt');
dg_ratio_nuc_right_file  = strcat(feature_dir,'/',prefix,'_',ids,'_R_DG_t2wt1wratio.txt');
sub_ratio_nuc_left_file  = strcat(feature_dir,'/',prefix,'_',ids,'_L_SUB_t2wt1wratio.txt');
sub_ratio_nuc_right_file = strcat(feature_dir,'/',prefix,'_',ids,'_R_SUB_t2wt1wratio.txt');

ca_ratio_nuc_left   = zeros(length(ids),10242);
ca_ratio_nuc_right  = zeros(length(ids),10242);
dg_ratio_nuc_left   = zeros(length(ids),5762);
dg_ratio_nuc_right  = zeros(length(ids),5762);
sub_ratio_nuc_left  = zeros(length(ids),5762);
sub_ratio_nuc_right = zeros(length(ids),5762);

for i=1:length(ids)

    disp(ids{i});
    ca_ratio_nuc_left(i,:)   = SurfStatReadData(ca_ratio_nuc_left_file{i});
    ca_ratio_nuc_right(i,:)  = SurfStatReadData(ca_ratio_nuc_right_file{i});
    dg_ratio_nuc_left(i,:)   = SurfStatReadData(dg_ratio_nuc_left_file{i});
    dg_ratio_nuc_right(i,:)  = SurfStatReadData(dg_ratio_nuc_right_file{i});
    sub_ratio_nuc_left(i,:)  = SurfStatReadData(sub_ratio_nuc_left_file{i});
    sub_ratio_nuc_right(i,:) = SurfStatReadData(sub_ratio_nuc_right_file{i});

end

%% z-score with respect to controls

% Global local volumes
columvol = [ca_columvol_left,sub_columvol_left,dg_columvol_left,ca_columvol_right,sub_columvol_right,dg_columvol_right];
columvol = SurfStatSmooth(columvol, template, 3);
[z_columvol_controls,   mu_columvol_controls,   std_columvol_controls]   = zscore(columvol(controls_group,:));
z_columvol_patients = (columvol(~controls_group,:) - repmat(mu_columvol_controls,[size(columvol(~controls_group,:),1),1]))./repmat(std_columvol_controls,[size(columvol(~controls_group,:),1),1]);

% T2 signal
t2_signal = [ca_t2_left,sub_t2_left,dg_t2_left,ca_t2_right,sub_t2_right,dg_t2_right];
t2_signal = SurfStatSmooth(t2_signal, template, 3);
[z_t2_signal_controls, mu_t2_signal_controls, std_t2_signal_controls] = zscore(t2_signal(controls_group,:));
z_t2_signal_patients = (t2_signal(~controls_group,:) - repmat(mu_t2_signal_controls,[size(t2_signal(~controls_group,:),1),1]))./repmat(std_t2_signal_controls,[size(t2_signal(~controls_group,:),1),1]);

% Ratio NUC Signal
ratio_nuc_signal = [ca_ratio_nuc_left,sub_ratio_nuc_left,dg_ratio_nuc_left,ca_ratio_nuc_right,sub_ratio_nuc_right,dg_ratio_nuc_right];
ratio_nuc_signal = SurfStatSmooth(ratio_nuc_signal, template, 3);
[z_ratio_nuc_signal_controls, mu_ratio_nuc_signal_controls, std_ratio_nuc_signal_controls] = zscore(ratio_nuc_signal(controls_group,:));
z_ratio_nuc_signal_patients = (ratio_nuc_signal(~controls_group,:) - repmat(mu_ratio_nuc_signal_controls,[size(ratio_nuc_signal(~controls_group,:),1),1]))./repmat(std_ratio_nuc_signal_controls,[size(ratio_nuc_signal(~controls_group,:),1),1]);

% Assymetry columnar volumes
columvol_assymetry = 2*(columvol(:,1:size(columvol,2)/2) - columvol(:,size(columvol,2)/2+1:end))./(columvol(:,1:size(columvol,2)/2) + columvol(:,size(columvol,2)/2+1:end));
[z_columvol_ass_controls, mu_columvol_ass_controls, std_columvol_ass_controls] = zscore(columvol_assymetry(controls_group,:));
z_columvol_ass_patients = (columvol_assymetry(~controls_group,:) - repmat(mu_columvol_ass_controls,[size(columvol_assymetry(~controls_group,:),1),1]))./repmat(std_columvol_ass_controls,[size(columvol_assymetry(~controls_group,:),1),1]);

% Assymetry T2 Signal
t2_signal_assymetry = 2*(t2_signal(:,1:size(t2_signal,2)/2) - t2_signal(:,size(t2_signal,2)/2+1:end))./(t2_signal(:,1:size(t2_signal,2)/2) + t2_signal(:,size(t2_signal,2)/2+1:end));
[z_t2_signal_ass_controls, mu_t2_signal_ass_controls, std_t2_signal_ass_controls] = zscore(t2_signal_assymetry(controls_group,:));
z_t2_signal_ass_patients = (t2_signal_assymetry(~controls_group,:) - repmat(mu_t2_signal_ass_controls,[size(t2_signal_assymetry(~controls_group,:),1),1]))./repmat(std_t2_signal_ass_controls,[size(t2_signal_assymetry(~controls_group,:),1),1]);

% Assymetry Ratio NUC Signal
ratio_nuc_signal_assymetry = 2*(ratio_nuc_signal(:,1:size(ratio_nuc_signal,2)/2) - ratio_nuc_signal(:,size(ratio_nuc_signal,2)/2+1:end))./(ratio_nuc_signal(:,1:size(ratio_nuc_signal,2)/2) + ratio_nuc_signal(:,size(ratio_nuc_signal,2)/2+1:end));
[z_ratio_nuc_signal_ass_controls, mu_ratio_nuc_signal_ass_controls, std_ratio_nuc_signal_ass_controls] = zscore(ratio_nuc_signal_assymetry(controls_group,:));
z_ratio_nuc_signal_ass_patients = (ratio_nuc_signal_assymetry(~controls_group,:) - repmat(mu_ratio_nuc_signal_ass_controls,[size(ratio_nuc_signal_assymetry(~controls_group,:),1),1]))./repmat(std_ratio_nuc_signal_ass_controls,[size(ratio_nuc_signal_assymetry(~controls_group,:),1),1]);

%% Switch from left right to ipsi contra

z_columvol_patients_ipsi_contra = zeros(size(z_columvol_patients));
z_t2_signal_patients_ipsi_contra = zeros(size(z_t2_signal_patients));
z_ratio_nuc_signal_patients_ipsi_contra = zeros(size(z_ratio_nuc_signal_patients));

z_columvol_ass_patients_ipsi_contra = zeros(size(z_columvol_ass_patients));
z_t2_signal_ass_patients_ipsi_contra = zeros(size(z_t2_signal_ass_patients));
z_ratio_nuc_signal_ass_patients_ipsi_contra = zeros(size(z_ratio_nuc_signal_ass_patients));

for i = 1:size(z_columvol_patients,1)

    if strcmp(lateralization_truth_training{i},'LTLE')
    
        z_columvol_patients_ipsi_contra(i,:) = z_columvol_patients(i,:);
        z_t2_signal_patients_ipsi_contra(i,:) = z_t2_signal_patients(i,:);
        z_ratio_nuc_signal_patients_ipsi_contra(i,:) = z_ratio_nuc_signal_patients(i,:);
        
        % z score of assymetries
        z_columvol_ass_patients_ipsi_contra(i,:) = z_columvol_ass_patients(i,:);
        z_t2_signal_ass_patients_ipsi_contra(i,:) = z_t2_signal_ass_patients(i,:);
        z_ratio_nuc_signal_ass_patients_ipsi_contra(i,:) = z_ratio_nuc_signal_ass_patients(i,:);
        
    else
        
        % z score of raw values
        z_columvol_patients_ipsi_contra(i,1:size(z_columvol_patients,2)/2) = z_columvol_patients(i,size(z_columvol_patients,2)/2+1:end);
        z_columvol_patients_ipsi_contra(i,size(z_columvol_patients,2)/2+1:end) = z_columvol_patients(i,1:size(z_columvol_patients,2)/2);
        
        z_t2_signal_patients_ipsi_contra(i,1:size(z_t2_signal_patients,2)/2) = z_t2_signal_patients(i,size(z_t2_signal_patients,2)/2+1:end);
        z_t2_signal_patients_ipsi_contra(i,size(z_t2_signal_patients,2)/2+1:end) = z_t2_signal_patients(i,1:size(z_t2_signal_patients,2)/2);
        
        z_ratio_nuc_signal_patients_ipsi_contra(i,1:size(z_ratio_nuc_signal_patients,2)/2) = z_ratio_nuc_signal_patients(i,size(z_ratio_nuc_signal_patients,2)/2+1:end);
        z_ratio_nuc_signal_patients_ipsi_contra(i,size(z_ratio_nuc_signal_patients,2)/2+1:end) = z_ratio_nuc_signal_patients(i,1:size(z_ratio_nuc_signal_patients,2)/2);
        
        % z score of assymetries
        z_columvol_ass_patients_ipsi_contra(i,:) = -z_columvol_ass_patients(i,:);
        z_t2_signal_ass_patients_ipsi_contra(i,:) = -z_t2_signal_ass_patients(i,:);
        z_ratio_nuc_signal_ass_patients_ipsi_contra(i,:) = -z_ratio_nuc_signal_ass_patients(i,:);
        
    end
    
end
