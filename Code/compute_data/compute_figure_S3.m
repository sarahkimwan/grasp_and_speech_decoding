%% compute data for figure S3 of the manuscript

%% Important: run code while being in folder 'grasp_and_speech_decoding'

clc
clear all 
close all

%% changeable parameters

dataFolder = [pwd '\Data\NeuralData\PCA_plot\'];
saveFolder = [pwd '\Data\CrossPhaseNeuronDropping\Speech'];

flagSaveData = false;
%if true  - saves computed data under saveFolder
%if false - does not save data

%%
testing_phases = {'Cue', 'Action'};

for i = 1:length(testing_phases)
[resultsPC{i},CueNamesPC{i}] = analysis.compute_NeuronDroppingCurve('SMG',...
    testing_phases{i},[dataFolder '\Table_sorting_aligned_thr_-4.5_SpokenColors.mat'],...
    {'Blue','Green','Yellow','Gray','Brown'});
end 

if flagSaveData
    save_sub_name = ['s2_SMG_PC_20_ColorSpeaking'];
    save_name = [saveFolder,save_sub_name];
    save(save_name, 'resultsPC', 'CueNamesPC');
end
clear resultsPC CueNamesPC

%SMG Grasp speaking
testing_phases = {'Cue', 'Action'};
for i = 1:length(testing_phases)

[resultsPC{i},CueNamesPC{i}] = analysis.compute_NeuronDroppingCurve('SMG',...
    testing_phases{i},[dataFolder '\Table_sorting_aligned_thr_-4.5_SpokenGrasps.mat'],...
    {'Lateral', 'WritingTripod', 'MediumWrap', 'PalmarPinch','Sphere3Finger'});
end 

if flagSaveData
    save_sub_name = ['s2_SMG_PC_20_GraspSpeaking'];
    save_name = [saveFolder,save_sub_name];
    save(save_name, 'resultsPC', 'CueNamesPC');
end 
clear resultsPC CueNamesPC

