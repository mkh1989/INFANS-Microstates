%% description
% This script reads the data and filters them into .
%
% NOTE: Some lines must be changed based on your data path.

%% initialization
close all; clear; clc;

subjectList   = setdiff(1:68, [10 19 24:25 45 59:60 66]); % excludes 8 neonates
dataPath      = '.\Data\';
savePath      = '.\Data\';

% new sampling ferquency
Fs    = 100;           % sampling rate
[p,q] = rat(Fs / 250); % for resampling

%% designing band pass filter (0.3-25 Hz)
lpButter = designfilt('lowpassiir', 'FilterOrder', 7, ...
'HalfPowerFrequency', 45, 'SampleRate', 250, 'DesignMethod', 'butter');     
hpButter = designfilt('highpassiir', 'FilterOrder', 7, ...
'HalfPowerFrequency', 0.15, 'SampleRate', 250, 'DesignMethod', 'butter');

%% band-pass filtering
for dataID = 1:length(subjectList)

    % loads .SET files for each subject 
    AS_Filename = ['AS_' num2str(subjectList(dataID)) '.set'];
    QS_Filename = ['QS_' num2str(subjectList(dataID)) '.set'];
    EEG_AS = pop_loadset([dataPath AS_Filename]);
    EEG_QS = pop_loadset([dataPath QS_Filename]);
    
    % changes the montage to the common average montage
    EEG_AS.data = infans_convert_to_averge_montage(EEG_AS.data);
    EEG_QS.data = infans_convert_to_averge_montage(EEG_QS.data);

    % removes the DC part from all EEG signals (based on Anton's suggestion)
    eeg_AS = eeg_AS - mean(eeg_AS, 2); % channels X samples
    eeg_QS = eeg_QS - mean(eeg_QS, 2); % channels X samples
    
    % filters, resmaples, and selects 3 minutes
    filtered_AS = filtfilt(hpButter,filtfilt(lpButter,eeg_AS')); % samples X channels
    filtered_QS = filtfilt(hpButter,filtfilt(lpButter,eeg_QS')); % samples X channels
    
    filtered_AS = resample(filtered_AS,p,q); % chnages Fs from 250 to 100 Hz
    filtered_QS = resample(filtered_QS,p,q); % chnages Fs from 250 to 100 Hz
    
    filtered_AS = filtered_AS(1: Fs * 180, :); % selects just 3 minutes
    filtered_QS = filtered_QS(1: Fs * 180, :); % selects just 3 minutes
    
    % save
    EEG = pop_saveset(EEG_AS,'filename',AS_Filename,'filepath',savePath,'savemode','onefile');
    EEG = pop_saveset(EEG_QS,'filename',QS_Filename,'filepath',savePath,'savemode','onefile');
     
    clear eeg_AS eeg_QS EEG
end

