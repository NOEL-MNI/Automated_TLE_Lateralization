%% Author: Benoit Caldairou, PhD
%% Mail: benoit.caldairou@mcgill.ca

%% Set up basic variables
n=100;
type = 'diagquadratic';
nFolds = 5;

%% Set up final ground truth and outputs
% Binarize lateralization
lateralization_bin_training = logical(strcmp(lateralization_truth_training,'LTLE'));

% Raw data labels
finalClassify_colvol_training       = zeros(length(lateralization_bin_training),n);
posterior_colvol_training           = zeros(length(lateralization_bin_training),n);

finalClassify_t2_training           = zeros(length(lateralization_bin_training),n);
posterior_t2_training              = zeros(length(lateralization_bin_training),n);

finalClassify_ratio_nuc_training    = zeros(length(lateralization_bin_training),n);
posterior_ratio_nuc_training        = zeros(length(lateralization_bin_training),n);

finalClassify_mult_t2_ratio_training    = zeros(length(lateralization_bin_training),n);
posterior_mult_t2_ratio_training        = zeros(length(lateralization_bin_training),n);

%% Cell arrays to store classifiers
lda_colvol          = cell(n,nFolds);
lda_t2              = cell(n,nFolds);
lda_ratio_nuc       = cell(n,nFolds);
lda_mult_t2_ratio   = cell(n,nFolds);

%% SPAM of ROIs
roi_spam_colvol_training = zeros(n,nFolds,21766);
roi_spam_t2_training = zeros(n,nFolds,21766);
roi_spam_ratio_nuc_training = zeros(n,nFolds,21766);
roi_spam_t2_ratio_training = zeros(n,nFolds,21766);

%% Run the validation
for i=1:n
    
    disp(i);
    partition = cvpartition(lateralization_bin_training,'KFold',nFolds);
    
    for k = 1:partition.NumTestSets
        
        disp(k);
        
        % ------------------------- ColVol ------------------------------ %
        % t-test to get t-maps beween ipsi and contra for training set
        [~,~,~,stats_nohs] = ttest2(z_columvol_patients_ipsi_contra(partition.training(k),1:size(z_columvol_patients,2)/2),z_columvol_patients_ipsi_contra(partition.training(k),size(z_columvol_patients,2)/2+1:end), 0.025, 'both','unequal');
        tmap = abs(stats_nohs.tstat);
        
        % Looks for optimal t in training set
        training_data  = z_columvol_ass_patients(partition.training(k),:);
        training_group = lateralization_bin_training(partition.training(k),:);
        t_optimal = optimal_thres_lookup(tmap,training_data,training_group);
        
        % Build data for this specific t
        roi = tmap >= t_optimal;
        roi_spam_colvol_training(i,k,:) = roi;
        database = mean(z_columvol_ass_patients(:,roi),2);
        
        % Build the model and classify
        lda_colvol{i,k} = fitcdiscr(database(partition.training(k)),lateralization_bin_training(partition.training(k)),'DiscrimType','diagquadratic');
        [finalClassify_colvol_training(partition.test(k),i),tmp_posterior] = predict(lda_colvol{i,k},database(partition.test(k)));
        posterior_colvol_training(partition.test(k),i) = tmp_posterior(:,2);
                       
        % ------------------------ T2 SIGNAL ---------------------------- %
        % t-test to get t-maps beween ipsi and contra for training set
        [~,~,~,stats_nohs] = ttest2(z_t2_signal_patients_ipsi_contra(partition.training(k),1:size(z_t2_signal_patients,2)/2),z_t2_signal_patients_ipsi_contra(partition.training(k),size(z_t2_signal_patients,2)/2+1:end), 0.025, 'both','unequal');
        tmap = abs(stats_nohs.tstat);
        
        % Looks for optimal t in training set
        training_data  = z_t2_signal_ass_patients(partition.training(k),:);
        training_group = lateralization_bin_training(partition.training(k),:);
        t_optimal = optimal_thres_lookup(tmap,training_data,training_group);
        
        % Build data for this specific t
        roi = tmap >= t_optimal;
        roi_spam_t2_training(i,k,:) = roi;
        database = mean(z_t2_signal_ass_patients(:,roi),2);
        
        % Build the model and classify
        lda_t2{i,k} = fitcdiscr(database(partition.training(k)),lateralization_bin_training(partition.training(k)),'DiscrimType','diagquadratic');
        [finalClassify_t2_training(partition.test(k),i),tmp_posterior] = predict(lda_t2{i,k},database(partition.test(k)));
        posterior_t2_training(partition.test(k),i) = tmp_posterior(:,2);
        
        % -------------------- RATIO NUC SIGNAL ------------------------- %
        % t-test to get t-maps beween ipsi and contra for training set
        [~,~,~,stats_nohs] = ttest2(z_ratio_nuc_signal_patients_ipsi_contra(partition.training(k),1:size(z_ratio_nuc_signal_patients,2)/2),z_ratio_nuc_signal_patients_ipsi_contra(partition.training(k),size(z_ratio_nuc_signal_patients,2)/2+1:end), 0.025, 'both','unequal');
        tmap = abs(stats_nohs.tstat);
        
        % Looks for optimal t in training set
        training_data  = z_ratio_nuc_signal_ass_patients(partition.training(k),:);
        training_group = lateralization_bin_training(partition.training(k),:);
        t_optimal = optimal_thres_lookup(tmap,training_data,training_group);
        
        % Build data for this specific t
        roi = tmap >= t_optimal;
        roi_spam_ratio_nuc_training(i,k,:) = roi;
        database = mean(z_ratio_nuc_signal_ass_patients(:,roi),2);
        
        % Build the model and classify
        lda_ratio_nuc{i,k} = fitcdiscr(database(partition.training(k)),lateralization_bin_training(partition.training(k)),'DiscrimType','diagquadratic');
        [finalClassify_ratio_nuc_training(partition.test(k),i),tmp_posterior] = predict(lda_ratio_nuc{i,k},database(partition.test(k)));
        posterior_ratio_nuc_training(partition.test(k),i) = tmp_posterior(:,2);
 
        % ------------- MULTIVARIATE T2 Ratio ----------------------- %
        clear forT2test_z_ipsi_contra_columnarvolume_training forT2test_z_ipsi_contra_ratiosignal_training forT2test_z_ipsi_contra_all_training;
        
        % t-test to get t-maps beween ipsi and contra for training set
        forT2test_z_ipsi_contra_columnarvolume_training = [z_columvol_patients_ipsi_contra(partition.training(k),1:size(z_columvol_patients,2)/2);z_columvol_patients_ipsi_contra(partition.training(k),size(z_columvol_patients,2)/2+1:end)];
        forT2test_z_ipsi_contra_t2_signal_training = [z_t2_signal_patients_ipsi_contra(partition.training(k),1:size(z_columvol_patients,2)/2);z_t2_signal_patients_ipsi_contra(partition.training(k),size(z_t2_signal_patients,2)/2+1:end)];
        forT2test_z_ipsi_contra_ratiosignalnuc_training = [z_ratio_nuc_signal_patients_ipsi_contra(partition.training(k),1:size(z_ratio_nuc_signal_patients,2)/2);z_ratio_nuc_signal_patients_ipsi_contra(partition.training(k),size(z_ratio_nuc_signal_patients,2)/2+1:end)];
        forT2test_z_ipsi_contra_all_training = cat(3,forT2test_z_ipsi_contra_t2_signal_training,forT2test_z_ipsi_contra_ratiosignalnuc_training);
        
        groups_ipsi_contra = {};
        groups_ipsi_contra(1:size(forT2test_z_ipsi_contra_columnarvolume_training,1)/2) = {'IPSI'};
        groups_ipsi_contra(size(forT2test_z_ipsi_contra_columnarvolume_training,1)/2+1:size(forT2test_z_ipsi_contra_columnarvolume_training,1)) = {'CONTRA'};
        GROUPS_IPSI_CONTRA = term(cellstr(groups_ipsi_contra));
        
        % Put a very simple model
        M = 1 + GROUPS_IPSI_CONTRA;
        % Hotelling T2 to get T-maps beween ipsi and contra
        slm_all = SurfStatLinMod(forT2test_z_ipsi_contra_all_training, M, template_uni); 
        slm_all = SurfStatT(slm_all,GROUPS_IPSI_CONTRA.IPSI-GROUPS_IPSI_CONTRA.CONTRA);% A few more variables
        
        tmap = abs(slm_all.t);
        
        clear training_data database;
        
        % Looks for optimal t in training set
        training_data(:,:,1)  = z_t2_signal_ass_patients(partition.training(k),:);
        training_data(:,:,2)  = z_ratio_nuc_signal_ass_patients(partition.training(k),:);
        training_group = lateralization_bin_training(partition.training(k),:);
        t_optimal = optimal_thres_lookup(tmap,training_data,training_group);
        
        % Build data for this specific t
        roi = tmap >= t_optimal;
        roi_spam_t2_ratio_training(i,k,:) = roi;
        database(:,1) = mean(z_t2_signal_ass_patients(:,roi),2);
        database(:,2) = mean(z_ratio_nuc_signal_ass_patients(:,roi),2);
       
        % Build the model and classify
        lda_mult_t2_ratio{i,k} = fitcdiscr(database(partition.training(k),:),lateralization_bin_training(partition.training(k)),'DiscrimType','linear');
        [finalClassify_mult_t2_ratio_training(partition.test(k),i),tmp_posterior] = predict(lda_mult_t2_ratio{i,k},database(partition.test(k),:));
        posterior_mult_t2_ratio_training(partition.test(k),i) = tmp_posterior(:,2);

    end

end

%% Check the results Surface Based
clear result_everyone_t2_training;
for i = 1:n

    result_everyone_t2_training(:,i) = finalClassify_t2_training(:,i) == lateralization_bin_training;
    
end

clear result_everyone_colvol_training;
for i = 1:n

    result_everyone_colvol_training(:,i) = finalClassify_colvol_training(:,i) == lateralization_bin_training;
    
end

clear result_everyone_ratio_nuc_training;
for i = 1:n

    result_everyone_ratio_nuc_training(:,i) = finalClassify_ratio_nuc_training(:,i) == lateralization_bin_training;
    
end

clear result_everyone_mult_t2_ratio_training;
for i = 1:n

    result_everyone_mult_t2_ratio_training(:,i) = finalClassify_mult_t2_ratio_training(:,i) == lateralization_bin_training;
    
end

%% Evaluate Global performance with surface based features
disp('------- All Patients ---------');

n_nonHS = length(lateralization_bin_training);

disp('Individual Performance')
["IDs", "Truth", "ColVol", "T2", "Ratio Nuc", "T2 Ratio";
ids_patients, lateralization_truth_training, num2cell(sum(result_everyone_colvol_training,2)), num2cell(sum(result_everyone_t2_training,2)), num2cell(sum(result_everyone_ratio_nuc_training,2)), num2cell(sum(result_everyone_mult_t2_ratio_training,2))]

disp('Average and Standard deviation')
["Structure", "Mean Performance", "Std Performance", "% Mean", "% Std";
"ColVol", num2cell(mean(sum(result_everyone_colvol_training))), num2cell(std(sum(result_everyone_colvol_training))), num2cell(mean(sum(result_everyone_colvol_training))*100/n_nonHS), num2cell(std(sum(result_everyone_colvol_training))*100/n_nonHS);
"T2", num2cell(mean(sum(result_everyone_t2_training))), num2cell(std(sum(result_everyone_t2_training))), num2cell(mean(sum(result_everyone_t2_training))*100/n_nonHS), num2cell(std(sum(result_everyone_t2_training))*100/n_nonHS);
"Ratio NUC", num2cell(mean(sum(result_everyone_ratio_nuc_training))), num2cell(std(sum(result_everyone_ratio_nuc_training))), num2cell(mean(sum(result_everyone_ratio_nuc_training))*100/n_nonHS), num2cell(std(sum(result_everyone_ratio_nuc_training))*100/n_nonHS);
"T2 Ratio", num2cell(mean(sum(result_everyone_mult_t2_ratio_training))), num2cell(std(sum(result_everyone_mult_t2_ratio_training))), num2cell(mean(sum(result_everyone_mult_t2_ratio_training))*100/n_nonHS), num2cell(std(sum(result_everyone_mult_t2_ratio_training))*100/n_nonHS);
]

disp('------- HS Part ---------');

n_nonHS = length(lateralization_bin_training(strcmp(hs_patients_training,'yes')));

disp('Individual Performance')
["IDs", "Truth", "ColVol", "T2", "Ratio Nuc", "T2 Ratio";
ids_patients(strcmp(hs_patients_training,'yes')), lateralization_truth_training(strcmp(hs_patients_training,'yes')), num2cell(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'yes'),:),2)), num2cell(sum(result_everyone_t2_training(strcmp(hs_patients_training,'yes'),:),2)), num2cell(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'yes'),:),2)), num2cell(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'yes'),:),2))]

disp('Average and Standard deviation')
["Structure", "Mean Performance", "Std Performance", "% Mean", "% Std";
"ColVol", num2cell(mean(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'yes'),:)))), num2cell(std(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'yes'),:)))), num2cell(mean(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'yes'),:)))*100/n_nonHS), num2cell(std(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'yes'),:)))*100/n_nonHS);
"T2", num2cell(mean(sum(result_everyone_t2_training(strcmp(hs_patients_training,'yes'),:)))), num2cell(std(sum(result_everyone_t2_training(strcmp(hs_patients_training,'yes'),:)))), num2cell(mean(sum(result_everyone_t2_training(strcmp(hs_patients_training,'yes'),:)))*100/n_nonHS), num2cell(std(sum(result_everyone_t2_training(strcmp(hs_patients_training,'yes'),:)))*100/n_nonHS);
"Ratio NUC", num2cell(mean(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'yes'),:)))), num2cell(std(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'yes'),:)))), num2cell(mean(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'yes'),:)))*100/n_nonHS), num2cell(std(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'yes'),:)))*100/n_nonHS);
"T2 Ratio", num2cell(mean(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'yes'),:)))), num2cell(std(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'yes'),:)))), num2cell(mean(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'yes'),:)))*100/n_nonHS), num2cell(std(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'yes'),:)))*100/n_nonHS);
]

disp('------- NonHS Part ---------');

n_nonHS = length(lateralization_bin_training(strcmp(hs_patients_training,'no')));

disp('Individual Performance')
["IDs", "Truth", "ColVol", "T2", "Ratio Nuc", "T2 Ratio";
ids_patients(strcmp(hs_patients_training,'no')), lateralization_truth_training(strcmp(hs_patients_training,'no')), num2cell(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'no'),:),2)), num2cell(sum(result_everyone_t2_training(strcmp(hs_patients_training,'no'),:),2)), num2cell(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'no'),:),2)), num2cell(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'no'),:),2))]

disp('Average and Standard deviation')
["Structure", "Mean Performance", "Std Performance", "% Mean", "% Std";
"ColVol", num2cell(mean(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'no'),:)))), num2cell(std(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'no'),:)))), num2cell(mean(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'no'),:)))*100/n_nonHS), num2cell(std(sum(result_everyone_colvol_training(strcmp(hs_patients_training,'no'),:)))*100/n_nonHS);
"T2", num2cell(mean(sum(result_everyone_t2_training(strcmp(hs_patients_training,'no'),:)))), num2cell(std(sum(result_everyone_t2_training(strcmp(hs_patients_training,'no'),:)))), num2cell(mean(sum(result_everyone_t2_training(strcmp(hs_patients_training,'no'),:)))*100/n_nonHS), num2cell(std(sum(result_everyone_t2_training(strcmp(hs_patients_training,'no'),:)))*100/n_nonHS);
"Ratio", num2cell(mean(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'no'),:)))), num2cell(std(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'no'),:)))), num2cell(mean(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'no'),:)))*100/n_nonHS), num2cell(std(sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'no'),:)))*100/n_nonHS);
"T2 Ratio", num2cell(mean(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'no'),:)))), num2cell(std(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'no'),:)))), num2cell(mean(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'no'),:)))*100/n_nonHS), num2cell(std(sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'no'),:)))*100/n_nonHS);
]

%% Significance of differences

disp('------- All Patients -------');
disp('Friedman Test');
all_results_training = [ sum(result_everyone_colvol_training)',... 
                sum(result_everyone_t2_training)',... 
                sum(result_everyone_ratio_nuc_training)',...
                sum(result_everyone_mult_t2_ratio_training)',...
                ];
            
[p_all_friedman_training,tbl,stats_all_friedman_training] = friedman(all_results_training);
c_all_friedman_training = multcompare(stats_all_friedman_training,'CType','bonferroni');

disp('------- HS Patients -------');
disp('Friedman Test');
hs_results_training = [  sum(result_everyone_colvol_training(strcmp(hs_patients_training,'yes'),:))',... 
                sum(result_everyone_t2_training(strcmp(hs_patients_training,'yes'),:))',... 
                sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'yes'),:))',...
                sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'yes'),:))',...
                ];
            
[p_hs_friedman_training,tbl,stats_hs_friedman_training] = friedman(hs_results_training);
c_hs_friedman_training = multcompare(stats_hs_friedman_training,'CType','bonferroni');

disp('------- Non HS Patients -------');
disp('Friedman Test');
nonhs_results_training = [ sum(result_everyone_colvol_training(strcmp(hs_patients_training,'no'),:))',... 
                  sum(result_everyone_t2_training(strcmp(hs_patients_training,'no'),:))',... 
                  sum(result_everyone_ratio_nuc_training(strcmp(hs_patients_training,'no'),:))',...
                  sum(result_everyone_mult_t2_ratio_training(strcmp(hs_patients_training,'no'),:))',...
                  ];
            
[p_nonhs_friedman_training,tbl,stats_nonhs_friedman_training] = friedman(nonhs_results_training);
c_nonhs_friedman = multcompare(stats_nonhs_friedman_training,'CType','bonferroni');

%% Confusion Matrices
matrix_colvol_training_HS = zeros(2,2,n);
matrix_colvol_training_nonHS = zeros(2,2,n);
matrix_t2_training_HS = zeros(2,2,n);
matrix_t2_training_nonHS = zeros(2,2,n);
matrix_ratio_nuc_training_HS = zeros(2,2,n);
matrix_ratio_nuc_training_nonHS = zeros(2,2,n);
matrix_t2_ratio_training_HS = zeros(2,2,n);
matrix_t2_ratio_training_nonHS = zeros(2,2,n);

for i = 1:n
    
    matrix_colvol_training_HS(:,:,i) = confusionmat(logical(finalClassify_colvol_training(strcmp(hs_patients_training,'yes'),i)),lateralization_bin_training(strcmp(hs_patients_training,'yes')),'Order',[1 0]);
    matrix_colvol_training_nonHS(:,:,i) = confusionmat(logical(finalClassify_colvol_training(strcmp(hs_patients_training,'no'),i)),lateralization_bin_training(strcmp(hs_patients_training,'no')),'Order',[1 0]);
    
    matrix_t2_training_HS(:,:,i) = confusionmat(logical(finalClassify_t2_training(strcmp(hs_patients_training,'yes'),i)),lateralization_bin_training(strcmp(hs_patients_training,'yes')),'Order',[1 0]);
    matrix_t2_training_nonHS(:,:,i) = confusionmat(logical(finalClassify_t2_training(strcmp(hs_patients_training,'no'),i)),lateralization_bin_training(strcmp(hs_patients_training,'no')),'Order',[1 0]);
    
    matrix_ratio_nuc_training_HS(:,:,i) = confusionmat(logical(finalClassify_ratio_nuc_training(strcmp(hs_patients_training,'yes'),i)),lateralization_bin_training(strcmp(hs_patients_training,'yes')),'Order',[1 0]);
    matrix_ratio_nuc_training_nonHS(:,:,i) = confusionmat(logical(finalClassify_ratio_nuc_training(strcmp(hs_patients_training,'no'),i)),lateralization_bin_training(strcmp(hs_patients_training,'no')),'Order',[1 0]);
    
    matrix_t2_ratio_training_HS(:,:,i) = confusionmat(logical(finalClassify_mult_t2_ratio_training(strcmp(hs_patients_training,'yes'),i)),lateralization_bin_training(strcmp(hs_patients_training,'yes')),'Order',[1 0]);
    matrix_t2_ratio_operated_nonHS(:,:,i) = confusionmat(logical(finalClassify_mult_t2_ratio_training(strcmp(hs_patients_training,'no'),i)),lateralization_bin_training(strcmp(hs_patients_training,'no')),'Order',[1 0]);
    

end

%% ROC curves
figure(1);

meanLineWidth = 4.0;
meanLineColor = '#ff0000';
meanLineStyle = '-';

individualLineColor = '#0080ff';
individualLineStyle = ':';

textFontSize = 24;

nb_ltle = sum(lateralization_bin_training);
nb_rtle = sum(~lateralization_bin_training);

nb_ltle_hs      = sum(lateralization_bin_training & strcmp(hs_patients_training,'yes'));
nb_ltle_nonhs   = sum(lateralization_bin_training & strcmp(hs_patients_training,'no'));

nb_rtle_hs      = sum(~lateralization_bin_training & strcmp(hs_patients_training,'yes'));
nb_rtle_nonhs   = sum(~lateralization_bin_training & strcmp(hs_patients_training,'no'));

% Colvol Subplot

clear -regex ltle_positive_* rtle_positive_* auc* ;

ax= subplot(2,2,1,'FontSize',textFontSize); hold on;
ax.XLabel.String = 'LTLE FPR'; ax.YLabel.String = 'LTLE TPR';

plot(0:0.001:1,0:0.001:1,'-k');

for repeat = 1:n

    tmp_posterior = posterior_colvol_training(:,repeat);
    t=0;
    
    for thres=0:0.001:1
    
        t = t+1;
        tmp_classify = tmp_posterior > thres;
        tmp_result = lateralization_bin_training == tmp_classify;
        
        tmp_result_hs = tmp_result(strcmp(hs_patients_training,'yes'));
        tmp_result_nonhs = tmp_result(strcmp(hs_patients_training,'no'));
        
        ltle_positive_colvol_training(repeat,t) = sum(tmp_result & lateralization_bin_training);
        rtle_positive_colvol_training(repeat,t) = sum(tmp_result & ~lateralization_bin_training);
        
        ltle_positive_colvol_training_hs(repeat,t) = sum(tmp_result_hs & lateralization_bin_training(strcmp(hs_patients_training,'yes')));
        rtle_positive_colvol_training_hs(repeat,t) = sum(tmp_result_hs & ~lateralization_bin_training(strcmp(hs_patients_training,'yes')));
        
        ltle_positive_colvol_training_nonhs(repeat,t) = sum(tmp_result_nonhs & lateralization_bin_training(strcmp(hs_patients_training,'no')));
        rtle_positive_colvol_training_nonhs(repeat,t) = sum(tmp_result_nonhs & ~lateralization_bin_training(strcmp(hs_patients_training,'no')));
        
        %ltle_positive_colvol_operated_hs(repeat,t) = ltle
    end

    ltle_positive_colvol_training(repeat,:) = ltle_positive_colvol_training(repeat,:)/nb_ltle;
    rtle_positive_colvol_training(repeat,:) = rtle_positive_colvol_training(repeat,:)/nb_rtle;
    
    ltle_positive_colvol_training_hs(repeat,:) = ltle_positive_colvol_training_hs(repeat,:)/nb_ltle_hs;
    rtle_positive_colvol_training_hs(repeat,:) = rtle_positive_colvol_training_hs(repeat,:)/nb_rtle_hs;
    
    ltle_positive_colvol_training_nonhs(repeat,:) = ltle_positive_colvol_training_nonhs(repeat,:)/nb_ltle_nonhs;
    rtle_positive_colvol_training_nonhs(repeat,:) = rtle_positive_colvol_training_nonhs(repeat,:)/nb_rtle_nonhs;
    
    auc_colvol_training(repeat)         = -trapz(1-rtle_positive_colvol_training(repeat,:),ltle_positive_colvol_training(repeat,:));
    auc_colvol_training_hs(repeat,:)    = -trapz(1-rtle_positive_colvol_training_hs(repeat,:),ltle_positive_colvol_training_hs(repeat,:));
    auc_colvol_training_nonhs(repeat,:) = -trapz(1-rtle_positive_colvol_training_nonhs(repeat,:),ltle_positive_colvol_training_nonhs(repeat,:));
    
    plot(1-rtle_positive_colvol_training_nonhs(repeat,:),ltle_positive_colvol_training_nonhs(repeat,:),'LineStyle',individualLineStyle,'Color', individualLineColor);

end

plot(mean(1-rtle_positive_colvol_training_nonhs),mean(ltle_positive_colvol_training_nonhs),'LineStyle',meanLineStyle,'Color',meanLineColor,'LineWidth',meanLineWidth);
text(0.35,0.1,['AUC = ', num2str(mean(auc_colvol_training_nonhs),'%0.2g')],'FontSize', textFontSize);

hold off;

% T2 Subplot

subplot(2,2,2,'FontSize',textFontSize); hold on;
plot(0:0.001:1,0:0.001:1,'-k');

for repeat = 1:n

    tmp_posterior = posterior_t2_training(:,repeat);
    t=0;
    
    for thres=0:0.001:1
    
        t = t+1;
        tmp_classify = tmp_posterior > thres;
        tmp_result = lateralization_bin_training == tmp_classify;
        
        tmp_result_hs = tmp_result(strcmp(hs_patients_training,'yes'));
        tmp_result_nonhs = tmp_result(strcmp(hs_patients_training,'no'));
        
        ltle_positive_t2_training(repeat,t) = sum(tmp_result & lateralization_bin_training);
        rtle_positive_t2_training(repeat,t) = sum(tmp_result & ~lateralization_bin_training);
        
        ltle_positive_t2_training_hs(repeat,t) = sum(tmp_result_hs & lateralization_bin_training(strcmp(hs_patients_training,'yes')));
        rtle_positive_t2_training_hs(repeat,t) = sum(tmp_result_hs & ~lateralization_bin_training(strcmp(hs_patients_training,'yes')));
        
        ltle_positive_t2_training_nonhs(repeat,t) = sum(tmp_result_nonhs & lateralization_bin_training(strcmp(hs_patients_training,'no')));
        rtle_positive_t2_training_nonhs(repeat,t) = sum(tmp_result_nonhs & ~lateralization_bin_training(strcmp(hs_patients_training,'no')));
        
    end

    ltle_positive_t2_training(repeat,:) = ltle_positive_t2_training(repeat,:)/nb_ltle;
    rtle_positive_t2_training(repeat,:) = rtle_positive_t2_training(repeat,:)/nb_rtle;
    
    ltle_positive_t2_training_hs(repeat,:) = ltle_positive_t2_training_hs(repeat,:)/nb_ltle_hs;
    rtle_positive_t2_training_hs(repeat,:) = rtle_positive_t2_training_hs(repeat,:)/nb_rtle_hs;
    
    ltle_positive_t2_training_nonhs(repeat,:) = ltle_positive_t2_training_nonhs(repeat,:)/nb_ltle_nonhs;
    rtle_positive_t2_training_nonhs(repeat,:) = rtle_positive_t2_training_nonhs(repeat,:)/nb_rtle_nonhs;
    
    auc_t2_training(repeat)         = -trapz(1-rtle_positive_t2_training(repeat,:),ltle_positive_t2_training(repeat,:));
    auc_t2_training_hs(repeat,:)    = -trapz(1-rtle_positive_t2_training_hs(repeat,:),ltle_positive_t2_training_hs(repeat,:));
    auc_t2_training_nonhs(repeat,:) = -trapz(1-rtle_positive_t2_training_nonhs(repeat,:),ltle_positive_t2_training_nonhs(repeat,:));
    
    plot(1-rtle_positive_t2_training_nonhs(repeat,:),ltle_positive_t2_training_nonhs(repeat,:),'LineStyle',individualLineStyle,'Color', individualLineColor);

end

plot(mean(1-rtle_positive_t2_training_nonhs),mean(ltle_positive_t2_training_nonhs),'LineStyle',meanLineStyle,'Color',meanLineColor,'LineWidth',meanLineWidth);
text(0.35,0.1,['AUC = ', num2str(mean(auc_t2_training_nonhs),'%0.2g')],'FontSize', textFontSize);

hold off;

% Ratio Subplot

subplot(2,2,3,'FontSize',textFontSize); hold on;
plot(0:0.001:1,0:0.001:1,'-k');

for repeat = 1:n

    tmp_posterior = posterior_ratio_nuc_training(:,repeat);
    t=0;
    
    for thres=0:0.001:1
    
        t = t+1;
        tmp_classify = tmp_posterior > thres;
        tmp_result = lateralization_bin_training == tmp_classify;
        
        tmp_result_hs = tmp_result(strcmp(hs_patients_training,'yes'));
        tmp_result_nonhs = tmp_result(strcmp(hs_patients_training,'no'));
        
        ltle_positive_ratio_nuc_training(repeat,t) = sum(tmp_result & lateralization_bin_training);
        rtle_positive_ratio_nuc_training(repeat,t) = sum(tmp_result & ~lateralization_bin_training);
        
        ltle_positive_ratio_nuc_training_hs(repeat,t) = sum(tmp_result_hs & lateralization_bin_training(strcmp(hs_patients_training,'yes')));
        rtle_positive_ratio_nuc_training_hs(repeat,t) = sum(tmp_result_hs & ~lateralization_bin_training(strcmp(hs_patients_training,'yes')));
        
        ltle_positive_ratio_nuc_training_nonhs(repeat,t) = sum(tmp_result_nonhs & lateralization_bin_training(strcmp(hs_patients_training,'no')));
        rtle_positive_ratio_nuc_training_nonhs(repeat,t) = sum(tmp_result_nonhs & ~lateralization_bin_training(strcmp(hs_patients_training,'no')));
        
        
    end

    ltle_positive_ratio_nuc_training(repeat,:) = ltle_positive_ratio_nuc_training(repeat,:)/nb_ltle;
    rtle_positive_ratio_nuc_training(repeat,:) = rtle_positive_ratio_nuc_training(repeat,:)/nb_rtle;
    
    ltle_positive_ratio_nuc_training_hs(repeat,:) = ltle_positive_ratio_nuc_training_hs(repeat,:)/nb_ltle_hs;
    rtle_positive_ratio_nuc_training_hs(repeat,:) = rtle_positive_ratio_nuc_training_hs(repeat,:)/nb_rtle_hs;
    
    ltle_positive_ratio_nuc_training_nonhs(repeat,:) = ltle_positive_ratio_nuc_training_nonhs(repeat,:)/nb_ltle_nonhs;
    rtle_positive_ratio_nuc_training_nonhs(repeat,:) = rtle_positive_ratio_nuc_training_nonhs(repeat,:)/nb_rtle_nonhs;
    
    auc_ratio_nuc_operated(repeat)          = -trapz(1-rtle_positive_ratio_nuc_training(repeat,:),ltle_positive_ratio_nuc_training(repeat,:));
    auc_ratio_nuc_operated_hs(repeat,:)     = -trapz(1-rtle_positive_ratio_nuc_training_hs(repeat,:),ltle_positive_ratio_nuc_training_hs(repeat,:));
    auc_ratio_nuc_operated_nonhs(repeat,:)  = -trapz(1-rtle_positive_ratio_nuc_training_nonhs(repeat,:),ltle_positive_ratio_nuc_training_nonhs(repeat,:));
    
    plot(1-rtle_positive_ratio_nuc_training_nonhs(repeat,:),ltle_positive_ratio_nuc_training_nonhs(repeat,:),'LineStyle',individualLineStyle,'Color', individualLineColor);

end

plot(mean(1-rtle_positive_ratio_nuc_training_nonhs),mean(ltle_positive_ratio_nuc_training_nonhs),'LineStyle',meanLineStyle,'Color',meanLineColor,'LineWidth',meanLineWidth);
text(0.35,0.1,['AUC = ', num2str(mean(auc_ratio_nuc_operated_nonhs),'%0.2g')],'FontSize', textFontSize);

hold off;

% T2/Ratio Subplot

subplot(2,2,4,'FontSize',textFontSize); hold on;
plot(0:0.001:1,0:0.001:1,'-k');

for repeat = 1:n

    tmp_posterior = posterior_mult_t2_ratio_training(:,repeat);
    t=0;
    
    for thres=0:0.001:1
    
        t = t+1;
        tmp_classify = tmp_posterior > thres;
        tmp_result = lateralization_bin_training == tmp_classify;
        
        tmp_result_hs = tmp_result(strcmp(hs_patients_training,'yes'));
        tmp_result_nonhs = tmp_result(strcmp(hs_patients_training,'no'));
        
        ltle_positive_t2_ratio_training(repeat,t) = sum(tmp_result & lateralization_bin_training);
        rtle_positive_t2_ratio_training(repeat,t) = sum(tmp_result & ~lateralization_bin_training);
        
        ltle_positive_t2_ratio_training_hs(repeat,t) = sum(tmp_result_hs & lateralization_bin_training(strcmp(hs_patients_training,'yes')));
        rtle_positive_t2_ratio_training_hs(repeat,t) = sum(tmp_result_hs & ~lateralization_bin_training(strcmp(hs_patients_training,'yes')));
        
        ltle_positive_t2_ratio_training_nonhs(repeat,t) = sum(tmp_result_nonhs & lateralization_bin_training(strcmp(hs_patients_training,'no')));
        rtle_positive_t2_ratio_training_nonhs(repeat,t) = sum(tmp_result_nonhs & ~lateralization_bin_training(strcmp(hs_patients_training,'no')));
        
    end

    ltle_positive_t2_ratio_training(repeat,:) = ltle_positive_t2_ratio_training(repeat,:)/nb_ltle;
    rtle_positive_t2_ratio_training(repeat,:) = rtle_positive_t2_ratio_training(repeat,:)/nb_rtle;
    
    ltle_positive_t2_ratio_training_hs(repeat,:) = ltle_positive_t2_ratio_training_hs(repeat,:)/nb_ltle_hs;
    rtle_positive_t2_ratio_training_hs(repeat,:) = rtle_positive_t2_ratio_training_hs(repeat,:)/nb_rtle_hs;
    
    ltle_positive_t2_ratio_training_nonhs(repeat,:) = ltle_positive_t2_ratio_training_nonhs(repeat,:)/nb_ltle_nonhs;
    rtle_positive_t2_ratio_training_nonhs(repeat,:) = rtle_positive_t2_ratio_training_nonhs(repeat,:)/nb_rtle_nonhs;
    
    auc_t2_ratio_training(repeat) = -trapz(1-rtle_positive_t2_ratio_training(repeat,:),ltle_positive_t2_ratio_training(repeat,:));
    auc_t2_ratio_training_hs(repeat,:)    = -trapz(1-rtle_positive_t2_ratio_training_hs(repeat,:),ltle_positive_t2_ratio_training_hs(repeat,:));
    auc_t2_ratio_training_nonhs(repeat,:) = -trapz(1-rtle_positive_t2_ratio_training_nonhs(repeat,:),ltle_positive_t2_ratio_training_nonhs(repeat,:));
    
    plot(1-rtle_positive_t2_ratio_training_nonhs(repeat,:),ltle_positive_t2_ratio_training_nonhs(repeat,:),'LineStyle',individualLineStyle,'Color', individualLineColor);

end

plot(mean(1-rtle_positive_t2_ratio_training_nonhs),mean(ltle_positive_t2_ratio_training_nonhs),'LineStyle',meanLineStyle,'Color',meanLineColor,'LineWidth',meanLineWidth);
text(0.35,0.1,['AUC = ', num2str(mean(auc_t2_ratio_training_nonhs),'%0.2g')],'FontSize', textFontSize);

hold off;

