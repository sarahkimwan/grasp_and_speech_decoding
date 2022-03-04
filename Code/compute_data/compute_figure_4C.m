%% compute data for cross-task classiifcation (figure 4E)

%% Important: run code while being in folder 'grasp_and_speech_decoding'

clc
clear all 
close all

%% changeable parameters

unit_region = 'SMG'; 
%Change between SMG, PMV and S1X
flagShuffleTest = false;
% if true, shuffle lables to compute shuffle distribution
flagSaveData = false; 
% save data
SavedData = [pwd, '\Data\NeuralData\SpeechDataset\']; %location of neural preprocessed data

%% 

number_repetition =1; 
if flagShuffleTest    
    number_repetition = 1000;
end 

task_names = {'MotorImagery','Grasps', 'Colors'};
classifier = 'diaglinear'; 
number_phases = 4; 

for n_task = 1:length(task_names)

    task_cue = task_names{n_task}; %Colors Grasps MotorImagery
    %load dataset
    cue_type = ['Table_sorting_aligned_thr_-4.5_Speaking_' task_cue '_ActionPhase_All_perTrials.mat'];
    data = load(fullfile(SavedData, cue_type));  
    Go_data = data.Go_data;
    Sessions = unique(Go_data.session_date);

    number_sessions = length(Sessions);
    labels = Go_data.GoLabels;

    %predeclare variables 
    k = 40;
    errTrain = ones(number_repetition,k,number_phases,number_sessions); 
    errTest = ones(number_repetition,k,number_phases,number_sessions);

    keep_session_idx_to_remove = [];
    flagPCA = true;
    for sessionIdx = 1:number_sessions
        disp(['Classifying session ' Sessions{sessionIdx} ]);
        idx_this_session = ismember(Go_data.session_date, Sessions(sessionIdx));

        if strcmp('SMG', unit_region)
            SessionData = Go_data.SMG_Go(idx_this_session,:);
        elseif strcmp('PMV', unit_region)
            SessionData = Go_data.PMV_Go(idx_this_session,:);
        elseif strcmp('S1X', unit_region)
            SessionData = Go_data.S1X_Go(idx_this_session,:);
        end

        %we remove session days that have no sorted units or less than 5 sorted
        %units 
        if nnz(cell2mat(cellfun(@(x) isempty(x), SessionData, 'UniformOutput', false))) >1 || unique(cell2mat(cellfun(@(x) size(x,2), SessionData, 'UniformOutput', false))) < 5
            keep_session_idx_to_remove = [keep_session_idx_to_remove, sessionIdx];
            warning(['Removing session ' Sessions{sessionIdx}])
            continue
        end

        labels_per_session = labels(idx_this_session);  
        time_phase_labels = Go_data.time_phase_labels(idx_this_session);

        for rep = 1:number_repetition 

            for n_phases = 1:number_phases
                %extract data per task phase
                data_per_phase = arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phases,:),1),SessionData,time_phase_labels, 'UniformOutput', false);
                %randomize labels to generate shuffle distribution
                if flagShuffleTest
                    labels_per_session = labels_per_session(randperm(length(labels_per_session)));
                end 

                labels_adapted = arrayfun(@(x,y) ones(size(x{1,1},1),1)*y ,data_per_phase,labels_per_session, 'UniformOutput', false);
                cv = cvpartition(length(data_per_phase), 'LeaveOut');

                for cvRun = 1:cv.NumTestSets %
                    %define training and testing data and labels 
                    trIdx = find(cv.training(cvRun));
                    teIdx = find(cv.test(cvRun));

                    data_training = cell2mat(data_per_phase(trIdx));
                    labels_train = cell2mat(labels_adapted(trIdx)); 
                    data_testing = cell2mat(data_per_phase(teIdx)); 
                    labels_test = cell2mat(labels_adapted(teIdx)); 

                    %perform pca
                    if flagPCA
                        [coeff, score, latent,~, explained] = pca(data_training); %The centered data can be reconstructed by SCORE*COEFF'.
                        %select PC's explaining 90% of variance
                        variance_c = cumsum(explained);
                        idx_90 = find(variance_c>90);
                        %project data onto PC's
                        PCA_data_train = data_training*coeff;
                        PCA_data_test = data_testing*coeff;
                        PCA_dimensions = idx_90(1); 

                        %if there are not many units, all are kept
                        if PCA_dimensions < 10
                            if idx_90(end) >= 10
                                PCA_dimensions = 10;
                            else
                                PCA_dimensions = idx_90(end);
                            end
                        end 

                        data_training = PCA_data_train(:,1:PCA_dimensions);    
                        data_testing = PCA_data_test(:,1:PCA_dimensions);
                        keepNbrPCAUnits(sessionIdx) = PCA_dimensions; 

                    end 
                    %train model 
                    model = fitcdiscr(data_training, labels_train, 'DiscrimType', classifier); 
                    %predict labels
                    predicted_training = predict(model, data_training);
                    predicted_testing = predict(model, data_testing);

                    errTrain(rep,cvRun,n_phases,sessionIdx) = analysis.classerror(labels_train, predicted_training); 
                    errTest(rep,cvRun,n_phases,sessionIdx) = analysis.classerror(labels_test, predicted_testing); 
                end 
                   clear LabelsRepRun;
                   clear LabelsPredRepRun;
            end 
        end 
    end
    % save data
    if flagSaveData
        ShuffleLabel = {'', 'Shuffle'};
        save([pwd '\Data\ClassificationComparisonMotorImagerySpeech\' unit_region '_' task_cue '_Error'   ShuffleLabel{flagShuffleTest+1}], ...
            'errTest')
    end 
end 
