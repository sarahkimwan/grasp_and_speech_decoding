%% Compute and plot PCR projections (figure S2 C)

%% Important: run code while being in folder 'grasp_and_speech_decoding'
clc
clear all
close all

%load motor imagery data
motor_imagery = load(fullfile(pwd, 'Data\NeuralData\PCA_plot','Table_sorting_aligned_thr_-4.5_MotorImagery.mat'));

unit_region = 'SMG';
data_phase = 4; %1 = ITI, 2 = Cue, 3 = Delay, 4 = Action 

if strcmp(unit_region, 'SMG')
    SessionsToInclude = {'20171003','20190404', '20190409',   '20190417',    '20190423',...
'20190509',    '20190510',    '20190522',    '20190528', '20190618', '20190911', '20190923', '20190924', '20190930','20191016'};
    nsp = 1;
elseif strcmp(unit_region, 'PMV')
    SessionsToInclude = {'20171003','20190404', '20190409',   '20190417', '20190522', '20190528', '20190618', ...
    '20190911', '20190923', '20190924', '20190930','20191016'};
    nsp = 2;
elseif strcmp(unit_region, 'S1X')
   SessionsToInclude = {'20190404',  '20190417',    '20190423',...
'20190509',    '20190510',    '20190522',    '20190528', '20190618', '20190911', '20190923', '20190924', '20190930','20191016'};
    nsp =3;
end 

%extract SMG Go trial data
go_data = motor_imagery.Go_data;
session_idx = logical((go_data.nsp == nsp).*ismember(go_data.session_date, SessionsToInclude));
%names of grasps
graspNames = go_data.Properties.VariableNames(5:9);

data_tmp = go_data(session_idx,5:9);
data = arrayfun(@(x) cell2mat(data_tmp{x,:}), 1:size(data_tmp,1), 'UniformOutput', false);
fig = figure('units','normalized','outerposition',[0 0 0.15 0.25]);

%% plot PCA projection of motor imagery (figure S2C panel 1)
for n_phase = data_phase 
    %extract data for corresponding task phase
    data_per_cue = cell2mat(arrayfun(@(x) mean(data{x}(go_data.time_phase_labels{x} == n_phase,:),1), 1:size(data,2), 'UniformOutput', false)');
    data_per_cue = data_per_cue';
    %zscore data
    data_per_cue = zscore(data_per_cue);
    %perform PCA
    [coeff, score, latent] = pca(data_per_cue);
    %get colors for each label
    colors_grasps = cell2mat(util.get_color_rgb_codes(graspNames)');
    %project data onto coefficients
    pca_data = data_per_cue*coeff; 

    Labels = reshape(repmat([1:5],8,1), [],1);
    %plot projection of first 2 principal components 
    gscatter(pca_data(:,1), pca_data(:,2),Labels, colors_grasps)
    xlabel('PCA1');
    ylabel('PCA2');
    %figure settings
    set(gca, 'fontsize', 15,'XTick',[],'Ytick', [] )
    title(['Motor Imagery N = 819'])
    xlim([ -15,15])
    ylim([-20,20])
end 

%% plot PCA projection of speaking data (figure S2C panel 2&3)

unit_region = 'SMG';
nsp = 1;
data_tasks = {'SpokenGrasps','SpokenColors',}; 

for n_task = 1:length(data_tasks)
    %extract speech data
    speech = load(fullfile(pwd, 'Data\NeuralData\PCA_plot',['Table_sorting_aligned_thr_-4.5_' data_tasks{n_task} '.mat']));
    go_data = speech.Go_data;

    session_idx = logical((go_data.nsp == nsp));

    %extract data for grasps or colors
    if strcmp(data_tasks{n_task}, 'SpokenColors')
        DataIdx = [11,13:16];
    else
        DataIdx = 5:9;
    end 
    %extract names
    graspNames = go_data.Properties.VariableNames(DataIdx);
    
    data_tmp = go_data(session_idx,DataIdx);
    data = arrayfun(@(x) cell2mat(data_tmp{x,:}), 1:size(data_tmp,1), 'UniformOutput', false);
    fig = figure('units','normalized','outerposition',[0 0 0.15 0.25]);

    for n_phase = data_phase
        %extract data for corresponding task phase
        data_per_cue = cell2mat(arrayfun(@(x) nanmean(data{x}(go_data.time_phase_labels{x} == n_phase,:),1), 1:size(data,2), 'UniformOutput', false)');
        data_per_cue = data_per_cue';
        %zscore data
        data_per_cue = zscore(data_per_cue);
        [coeff, score, latent] = pca(data_per_cue);
        %extract color
        colors_grasps = cell2mat(util.get_color_rgb_codes(graspNames)');
        pca_data = data_per_cue*coeff; 

        Labels = reshape(repmat([1:5],8,1), [],1);
        %plot projection of first 2 principal components 
        gscatter(pca_data(:,1), pca_data(:,2),Labels, colors_grasps)
        %figure settings
         set(gca, 'fontsize', 15,'XTick',[],'Ytick', [] )
        title([data_tasks{n_task} '  N = 252'])
        xlabel('PCA1');
        ylabel('PCA2');
    end 
end 


