function [tuned_combined_units, tuned_units_per_phase, tuned_units_per_bin, sum_phase, sum_bin,nbr_tuned_units_per_grasp,nbr_tuned_units_per_grasp_bin,p_per_phase, p_per_phase_orig] = get_tuned_units(Data,Labels, TimePhaseLabels, flag50msTuning)


%Data | size[number_trials,1]. Each cell contains a matrix of firing rate in 50ms time bins of the task x number of units
%Labels | size[number_of_trials]. Contains label for each trial
%TimePhaseLabels | labels of phase in 50ms time bins. 
%flag50msTuning | true or false | defines if tuning in 50ms time bins is computed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graspNames = {'Lateral','WritingTripod','MediumWrap','PalmarPinch','Sphere3Finger'};

%variables
alpha = 0.05; %defines tuning
number_bins = size(Data{1},1);
number_units = size(Data{1},2);
number_grasps = length(graspNames);
unique_phase_labels = unique(TimePhaseLabels{1});
number_phases = unique_phase_labels(end); %number of phases in task
fixed_trial_number = unique(histc(Labels, unique(Labels))); % the number of trials should be the same for all 

condition_names = ['ITI',graspNames];
condition = repmat(condition_names, [fixed_trial_number,1]);
condition = reshape(condition, [fixed_trial_number*(number_grasps+1),1]);

%presence vector X
for i = 1:number_grasps
    X(:,i) = ismember(condition, graspNames{i});
end

%predeclare variables
p_per_phase = ones(number_phases,number_grasps, number_units);
p_per_phase_orig = ones(number_phases, number_grasps,number_units);
p_per_bin = ones(number_bins, number_grasps, number_units);

ITI_DataAll =  cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== 1,:),1),Data,TimePhaseLabels, 'UniformOutput', false));

for n_phase = 1:number_phases
    %extract firing rate for each phase and average results
    data_per_phase = cell2mat(arrayfun(@(x,y) nanmean(x{1,1}(y{:}== n_phase,:),1),Data,TimePhaseLabels, 'UniformOutput', false));
    %loop through each unit
    for n_unit = 1:number_units

      data_per_trial = data_per_phase(:,n_unit);
      %order the trials per label
      data_per_trial_ordered = cell2mat(arrayfun(@(x) data_per_trial(Labels == x), unique(Labels), 'UniformOutput', false));
      %compute average ITI as baseline
      ITI_data = ones(fixed_trial_number,1)*mean(mean(ITI_DataAll(:,n_unit)));
      FR = [ITI_data; data_per_trial_ordered]; 
      %fit linear regression model
      mdl = fitlm(X,FR);
      %p values of model per grasp 
      p_per_phase(n_phase,:,n_unit) = mdl.Coefficients.pValue(2:end); 
      p_per_phase_orig(n_phase,:,n_unit) =  p_per_phase(n_phase,:,n_unit);

    end       
end 

%if true, compute tuning for each 50ms time bin 
if flag50msTuning 
    for n_bin = 1:number_bins
        data_per_bin = cell2mat(arrayfun(@(x,y) Data{x,1}(n_bin,:),1:size(Data,1), 'UniformOutput', false)');
        %loop through number of units 
        for n_unit = 1:number_units
            data_per_bin_per_trial = data_per_bin(:,n_unit);
            data_per_bin_ordered = cell2mat(arrayfun(@(x) data_per_bin_per_trial(Labels == x), unique(Labels), 'UniformOutput', false));
            %compute average ITI as baseline
            ITI_data = ones(fixed_trial_number,1)*mean(mean(ITI_DataAll(:,n_unit)));
            %fit linear regression model
            FR = [ITI_data; data_per_bin_ordered]; 
            mdl = fitlm(X,FR);
            p_per_bin(n_bin,:,n_unit) = mdl.Coefficients.pValue(2:end); 
        end
    end 
end 

%output variables
tuned_units_per_phase_per_grasp = p_per_phase < alpha;
nbr_tuned_units_per_grasp = sum(tuned_units_per_phase_per_grasp,3);

tuned_units_per_phase = logical(squeeze(sum(tuned_units_per_phase_per_grasp,2)))';
sum_phase = sum(tuned_units_per_phase);

tuned_units_per_bin = p_per_bin < 0.05;
nbr_tuned_units_per_grasp_bin = sum(tuned_units_per_bin,3);

sum_bin = sum(logical(squeeze(sum(tuned_units_per_bin,2))),2); 
tuned_combined_units = logical(sum(tuned_units_per_phase,2));    

end

