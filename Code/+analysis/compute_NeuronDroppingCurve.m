% Cross phase classification + neuron dropping curve 

function [result,cue_names] = compute_NeuronDroppingCurve(unit_region,TrainingPhase,TableName, GroupNames)

%classification parameters
k = 8; %cross-validation number
classifier = 'diaglinear'; %classifier type
flagPCA = true; % compute PCA
n_PC = 20; % number of PC's included in anaysis
num_of_iterations = 100; %number of iterations of neuron dropping curve

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

phase_cues = {'Cue', 'Action'};
cue_count  = length(phase_cues);
n = 0;
data_per_phase_all = cell(cue_count,1);

for testing_phase  = [2 4]
    n = n +1;
    %load data
    load(TableName)

    %extract data
    if strcmp(unit_region, 'all')
        region_idx = ones(size(Go_data.nsp));
    elseif strcmp(unit_region, 'AIP')|| strcmp(unit_region, 'SMG') || strcmp(unit_region, 'M1')
       region_idx = (Go_data.nsp == 1); 
    elseif strcmp(unit_region, 'BA5') || strcmp(unit_region, 'PMV') || strcmp(unit_region, 'PPC')
       region_idx = (Go_data.nsp == 2); 
    elseif strcmp(unit_region, 'S1X') 
        region_idx = (Go_data.nsp == 3);    
    end
    
    %include sessions
    day_idx = ismember(Go_data.session_date, SessionsToInclude);
    region_idx = logical(region_idx.*day_idx); 
   
    data_all = table2cell(Go_data(region_idx,ismember(Go_data.Properties.VariableNames, GroupNames)));

    time_phase_labels = Go_data.time_phase_labels{1};

    if testing_phase == 2 %extract cue phase data
        phase_idx_to_test = time_phase_labels == testing_phase;
        
        if contains(TableName, 'SpokenColors')
            time_trial = 0.05*(1:length(phase_idx_to_test));
            startIdx = find(time_trial>4,1);
            %shorten cue phase data when SpokenColors data is used
            disp('shorter cue phase for spoken color')
            phase_idx_to_test(startIdx:end) = 0;
        end 
    else %extract action phase data
        phase_idx_to_test = time_phase_labels == testing_phase;
    end 
    %transform data into 40 trials x number of features matrix 
    data_phase = cell2mat(cellfun(@(x) nanmean(x(phase_idx_to_test,:)), data_all,'UniformOutput', false))';
    data_per_phase_all{n} = data_phase;
end


%%%% CHOOSE HERE %%%%%

if strcmp(TrainingPhase, 'Cue')
    data1_phase = data_per_phase_all{1}; %Cue Phase Data
elseif strcmp(TrainingPhase, 'Action')
    data1_phase = data_per_phase_all{2}; %Action Phase Data
end

if isequal(data1_phase,data_per_phase_all{1})

    data1_name = 'Cue';
    data2_phase = data_per_phase_all{2};    
    data2_name = 'Action';

elseif isequal(data1_phase,data_per_phase_all{2}) 

    data1_name = 'Action';
    data2_phase = data_per_phase_all{1};
    data2_name = 'Cue';
end

%extract labels
labels_number = cell2mat(cellfun(@(x) util.image2class_simple(x), GroupNames, 'UniformOutput', false));
labels_all = reshape(repmat(labels_number,8,1), [],1);

%number of available units
max_feature_number = size(data_phase,2);
%save results in output variable
result = cell(1,num_of_iterations);

tic
%loop through number of iterations
parfor rep = 1:num_of_iterations

    err = nan(k,cue_count,  max_feature_number);
    disp(['Repetition ' num2str(rep) '/' num2str(num_of_iterations)])
    for num_features = 1:max_feature_number
        
        %completely randomize the available features each time
        features_rand = randperm(size(data1_phase,2));
        %keep 1 feature in loop 1, 2 features in loop 2, etc. 
        features_idx = features_rand(1:num_features);
        %select the data with randomly selected units
        data_cv = data1_phase(:,features_idx);
        data_2 = data2_phase(:,features_idx);
        %cross-validation partition
        cv = cvpartition(size(data1_phase,1), 'KFold', k);
        cv = repartition(cv); 
        %iterate over different data partitions
        for runNbr = 1:cv.NumTestSets 
            %define training and testing data and labels 
            trIdx = find(cv.training(runNbr));
            teIdx = find(cv.test(runNbr));
            
            labels_train = labels_all(trIdx);
            labels_test = labels_all(teIdx);
            
            DataTraining = data_cv(trIdx,:);
            DataTesting1 = data_cv(teIdx,:);
            DataTesting2 = data_2(teIdx,:);
            
            %perform PCA
            if flagPCA
                     %compute coefficient on training data
                     [coeff] = pca(DataTraining);                      
                      n_PC_max = n_PC;
                      
                      %If the number of PC is bigger than the number of
                      %features, we keep all the features 
                      if n_PC_max < size(DataTraining,2) && n_PC_max < size(coeff,2)
                        %apply coefficients on train and testing data
                        DataTraining = DataTraining*coeff(:,1:n_PC_max);
                        DataTesting1 = DataTesting1*coeff(:,1:n_PC_max);
                        DataTesting2 = DataTesting2*coeff(:,1:n_PC_max);
                      else
                        DataTraining = DataTraining*coeff;
                        DataTesting1 = DataTesting1*coeff;
                        DataTesting2 = DataTesting2*coeff;
                      end 
            end 
            %train classification model
            model = fitcdiscr(DataTraining, labels_train, 'DiscrimType', classifier);
            %predict and calculate classification error
            err(runNbr,1 ,num_features) = analysis.classerror(labels_test, predict(model, DataTesting1));
            err(runNbr,2 ,num_features) = analysis.classerror(labels_test, predict(model, DataTesting2));
 
        end 
        
    end 
    result{rep} = err;
end 

toc

cue_names = {data1_name, data2_name};

end


