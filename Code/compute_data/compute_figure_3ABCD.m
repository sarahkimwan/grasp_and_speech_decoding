%% compute data for figure 3 of the manuscript

%% Important: run code while being in folder 'grasp_and_speech_decoding'
addpath(genpath(pwd)); %add folder to search path 

%%
clc
clear all 
close all

%% changeable parameters

dataFolder = [pwd '\Data\NeuralData\PCA_plot\'];
saveFolder = [pwd '\Data\CrossPhaseNeuronDropping\MotorImagery'];

flagSaveData = false;
%if true  - saves computed data under saveFolder
%if false - does not save data
%% SMG
testing_phases = {'Cue', 'Action'};
for i = 1:length(testing_phases)
[resultsPC{i},CueNamesPC{i}] = analysis.compute_NeuronDroppingCurve('SMG',...
    testing_phases{i},[dataFolder '\Table_sorting_aligned_thr_-4.5_MotorImagery.mat'],...
    {'Lateral', 'WritingTripod', 'MediumWrap', 'PalmarPinch','Sphere3Finger'});
end 

if flagSaveData
    saveSubName = ['s2_SMG_PC_20'];
    save_name = [saveFolder,saveSubName];
    save(save_name, 'resultsPC', 'CueNamesPC');
end 
clear resultsPC CueNamesPC

%% PMV
for i = 1:length(testing_phases)
[resultsPC{i},CueNamesPC{i}] = analysis.compute_NeuronDroppingCurve('PMV',...
    testing_phases{i},[dataFolder '\Table_sorting_aligned_thr_-4.5_MotorImagery.mat'],...
    {'Lateral', 'WritingTripod', 'MediumWrap', 'PalmarPinch','Sphere3Finger'});
end 

if flagSaveData
    saveSubName = ['s2_PMV_PC_20'];
    save_name = [saveFolder,saveSubName];
    save(save_name, 'resultsPC', 'CueNamesPC');
end 
clear resultsPC CueNamesPC

%% S1
for i = 1:length(testing_phases)
[resultsPC{i},CueNamesPC{i}] = analysis.compute_NeuronDroppingCurve('S1X',...
    testing_phases{i},[dataFolder '\Table_sorting_aligned_thr_-4.5_MotorImagery.mat'],...
    {'Lateral', 'WritingTripod', 'MediumWrap', 'PalmarPinch','Sphere3Finger'});
end 

if flagSaveData
    saveSubName = ['s2_S1_PC_20'];
    save_name = [saveFolder,saveSubName];
    save(save_name, 'resultsPC', 'CueNamesPC');
end 
clear resultsPC CueNamesPC
