function [result,cue_names] = compute_NeuronDroppingCurve(subject,unit_region,TrainingPhase, varargin)

%clc
%clear all

[varargin,num_of_iterations] = Utilities.argkeyval('NbrIterations',varargin, 1000);
[varargin,k] = Utilities.argkeyval('number_of_folds',varargin, 4);
[varargin,classifier] = Utilities.argkeyval('classifier',varargin, 'diaglinear');
[varargin,flagPCA] = Utilities.argkeyval('PCA', varargin, 'true'); 
[varargin,NbrPC] = Utilities.argkeyval('NbrPC', varargin, 30); %Change the value of PC components 


Utilities.argempty(varargin);

PhaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
%TestingPhase = find(ismember(PhaseNames, TestingPhaseName));
%subject = 'n1';

%unit_region = 'PPC'; 
%num_of_iterations = 2;

if strcmp(unit_region, 'SMG')
    SessionsToInclude = {'20171003','20190404', '20190409',   '20190417',    '20190423',...
'20190509',    '20190510',    '20190522',    '20190528', '20190618', '20190911', '20190923', '20190924', '20190930','20191016'};
elseif strcmp(unit_region, 'PMV')
    SessionsToInclude = {'20171003','20190404', '20190409',   '20190417', '20190522', '20190528', '20190618', ...
    '20190911', '20190923', '20190924', '20190930','20191016'};
elseif strcmp(unit_region, 'S1X')
   SessionsToInclude = {'20190404',  '20190417',    '20190423',...
'20190509',    '20190510',    '20190522',    '20190528', '20190618', '20190911', '20190923', '20190924', '20190930','20191016'};

end



cues = {'Cue', 'Action'};
cue_count  = length(cues);
%TestingPhase = 2; %find(ismember(data.task.phaseNames, testing_phase_name));
n = 0;
for TestingPhase  = [2 4]
    n = n +1;
    
    %load(['D:\Users\Sarah\Documents\Saved_Data\Table_Data\' subject '_2019_aligned_correctly\Table_sorting_aligned_thr_-4.5_' cueType 'CuePerModality_2s.mat'])
    warning('MAJE SURE THIS IS CORRECT !!! !')
    load(['C:\Users\Sarah\Documents\Saved_Data\Table_Data\' subject '_2020_aligned_correctly\Table_sorting_aligned_thr_-4.5_ImageCue_all.mat'])

    
    grasps_to_test = {'Lateral', 'WritingTripod', 'MediumWrap', 'PalmarPinch','Sphere3Finger'}; 

    classifier = 'diaglinear';

    %logical array for unit region
    if strcmp(unit_region, 'all')
        region_idx = ones(size(Go_data.nsp));

    elseif strcmp(unit_region, 'AIP')|| strcmp(unit_region, 'SMG') || strcmp(unit_region, 'M1')
       region_idx = (Go_data.nsp == 1); 
    elseif strcmp(unit_region, 'BA5') || strcmp(unit_region, 'PMV') || strcmp(unit_region, 'PPC')
       region_idx = (Go_data.nsp == 2); 
    elseif strcmp(unit_region, 'S1X') 
        region_idx = (Go_data.nsp == 3);    
    end
    
    day_idx = ismember(Go_data.session_date, SessionsToInclude);

    region_idx = logical(region_idx.*day_idx); 
    
    data_all = table2cell(Go_data(region_idx,ismember(Go_data.Properties.VariableNames, grasps_to_test)));

    time_phase_labels = Go_data.time_phase_labels{1};

    %transform data into 40 trials x NBR features matrix 
    data_phase = cell2mat(cellfun(@(x) mean(x(time_phase_labels == TestingPhase,:)), data_all,'UniformOutput', false))';
    data_phase_all{n} = data_phase;
end


%%%% CHOOSE HERE %%%%%

if strcmp(TrainingPhase, 'Cue')
    data_phase1 = data_phase_all{1}; %Cue Phase Data
elseif strcmp(TrainingPhase, 'Action')
    data_phase1 = data_phase_all{2}; %Action Phase Data
end


% define phase 1 and 2 data and name
if isequal(data_phase1,data_phase_all{1})
    data1_name = 'Cue';
    data_phase2 = data_phase_all{2};    
    data2_name = 'Action';
elseif isequal(data_phase1,data_phase_all{2}) 
    data1_name = 'Action';
    data_phase2 = data_phase_all{1};
    data2_name = 'Cue';
end 


labels_number = cell2mat(cellfun(@(x) preproc.image2class_simple(x), grasps_to_test, 'UniformOutput', false));
labels_all = reshape(repmat(labels_number,8,1), [],1);

max_feature_number = size(data_phase,2);
result = cell(1,num_of_iterations);
k = 8;

tic
parfor rep = 1:num_of_iterations

%for rep = 1:num_of_iterations %for debugging
    %err = nan(k,num_phases,max_feature_number);
    err = nan(k,cue_count,  max_feature_number);
    disp(rep)
    for num_features = 1:max_feature_number
        
        %completely randomize the available features each time
        features_rand = randperm(size(data_phase1,2));
        %keep 1 feature in loop 1, 2 features in loop 2, etc. 
        features_idx = features_rand(1:num_features);
        %select the data with the selected features (= units)
        data_cv = data_phase1(:,features_idx);
        data_2 = data_phase2(:,features_idx);
          
        cv = cvpartition(size(data_phase1,1), 'KFold', k);
        cv = repartition(cv); %pas forcement utile
        
        for runNbr = 1:cv.NumTestSets %iterate over different data partitions
            %find training and testing data and labels indexes
            trIdx = find(cv.training(runNbr));
            teIdx = find(cv.test(runNbr));
            
            labels_train = labels_all(trIdx);
            
            labels_test = labels_all(teIdx);
            
            data_training = data_cv(trIdx,:);
            
            data_testing1 = data_cv(teIdx,:);
            data_testing2 = data_2(teIdx,:);
            
            if flagPCA
                     [coeff, ~, ~,~, ~] = pca(data_training);
                       idx_90 = NbrPC; %find(variance_c>90,1);
                      if idx_90 < size(data_training,2) && idx_90 < size(coeff,2)
                        data_training = data_training*coeff(:,1:idx_90);
                        data_testing1 = data_testing1*coeff(:,1:idx_90);
                        data_testing2 = data_testing2*coeff(:,1:idx_90);
                      else
                        data_training = data_training*coeff;
                        data_testing1 = data_testing1*coeff;
                        data_testing2 = data_testing2*coeff;
                      end 
                      %KeepPC = [KeepPC, idx_90];
            end 
            
            model = fitcdiscr(data_training, labels_train, 'DiscrimType', classifier);
            %size(DataTraining)

            %predicted_labels1 = predict(model, DataTesting1); 
            err(runNbr,1 ,num_features) = classification.classerror(labels_test, predict(model, data_testing1));
            err(runNbr,2 ,num_features) = classification.classerror(labels_test, predict(model, data_testing2));
 
        end 
        
    end 
    result{rep} = err;

end 

toc

err1 = cell2mat(cellfun(@(x) squeeze(x(:,1,:)),result', 'UniformOutput', false));
err2 = cell2mat(cellfun(@(x) squeeze(x(:,2,:)),result', 'UniformOutput', false));

cue_names = {data1_name, data2_name};

color_modalities = utile.get_color_rgb_codes({data1_name, data2_name});
figure(); 

plot(1-mean(err1), 'Color', color_modalities{1});
hold on 
plot(1-mean(err2),'Color', color_modalities{2});
hold on 
title([subject  ' - ' unit_region ' - Training with ' data1_name ' Phase'] )
disp('here')

end


