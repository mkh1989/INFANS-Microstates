%% Initialization
clear; close all; clc;

ExcludedSubjects = [10 19 24 25 45 59 60 66];
SubjectList      = setdiff(1:68,ExcludedSubjects); % excludes 8 neonates

% Defining Paths
AS_Group_Dir = '.\Active Sleep';
TA_Group_Dir = '.\Quiet Sleep';
SavePath     = '.\Results';

% Defining clustering parameters
ClustPars.MinClasses     = 3;         % minimum number of microstates
ClustPars.MaxClasses     = 15;        % maximum number of microstates
ClustPars.GFPPeaks       = 1;         % use GFP peaks as inputs of clustering algorithm
ClustPars.IgnorePolarity = 1;
ClustPars.MaxMaps        = inf;       % number of time sample data feeded to clusterig algorithm 
ClustPars.Restarts       = 20;        % number of repeatition for modified k-means 
ClustPars.UseAAHC        = false;     % if true, TAAGC will be applied, else modified k-menas.

% Defining fitting parameters
FitPars.b        = 30;  % the window size, if labels are smoothed
FitPars.lambda   = 1;   % the non-smoothness penalty factor
FitPars.PeakFit  = 1;   % true if fitting is only done on GFP peaks and the assignment is interpolated in between, false otherwise
FitPars.BControl = 1;

eeglab

%% Reading the data
AS_Group_Index = [];
TA_Group_Index = [];

% Read the data from AS group
for i = 1:length(SubjectList)
    EEG = pop_loadset([AS_Group_Dir '\AS_' num2str(SubjectList(i)) '.set']); % read file
    if size(EEG.data,2) > 30000                                              % reduce the length of data for computational costs
        EEG.data  = EEG.data(:,1:30000);
    end
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'gui','off');  % make this a new set
    ALLEEG                    = eeg_store(ALLEEG, EEG, CURRENTSET);       % store the thing
    AS_Group_Index            = [AS_Group_Index CURRENTSET];              % keep track of the groups
end

% Read the data from TA group
for i = 1:length(SubjectList)
    EEG = pop_loadset([TA_Group_Dir '\TA_' num2str(SubjectList(i)) '.set']); % read file
    if size(EEG.data,2) > 30000                                              % reduce the length of data for computational costs
        EEG.data = EEG.data(:,1:30000);
    end
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'gui','off');  % make this a new set
    ALLEEG                    = eeg_store(ALLEEG, EEG, CURRENTSET);       % store the thing
    TA_Group_Index            = [TA_Group_Index CURRENTSET];              % keep track of the groups
end

AllSubjects = [AS_Group_Index TA_Group_Index];

eeglab redraw

%% Microstates for each subject
% loop across all subjects to identify the individual clusters
for i = 1:numel(AllSubjects) 
    tmpEEG = eeg_retrieve(ALLEEG,AllSubjects(i));                                     % the EEG we want to work with
    fprintf(1,'Clustering dataset %s (%i/%i)\n',tmpEEG.setname,i,numel(AllSubjects)); % some info for the impatient user
    tmpEEG = pop_FindMSTemplates(tmpEEG, ClustPars);                                  % this is the actual clustering within subjects
    ALLEEG = eeg_store(ALLEEG, tmpEEG, AllSubjects(i));                               % done, we just need to store this
end

eeglab redraw

%% Combining the microstate maps across subjects
% The mean of AS group
EEG = pop_CombMSTemplates(ALLEEG, AS_Group_Index, 0, 0, 'GrandMean AS Group');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1, 'gui', 'off'); % make a new set
[ALLEEG,EEG, CURRENTSET]  = eeg_store(ALLEEG, EEG, CURRENTSET);                     % store it
Grand_Mean_AS_G_Index     = CURRENTSET;                                             % keep track of it

% The mean of TA group
EEG = pop_CombMSTemplates(ALLEEG, TA_Group_Index, 0, 0, 'GrandMean TA Group');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1, 'gui', 'off'); % make a new set
[ALLEEG,EEG, CURRENTSET]  = eeg_store(ALLEEG, EEG, CURRENTSET);                     % store it
Grand_Mean_TA_G_Index     = CURRENTSET;                                             % keep track of it

% The grand-grand mean, based on the two group means
EEG = pop_CombMSTemplates(ALLEEG, [Grand_Mean_AS_G_Index Grand_Mean_TA_G_Index], 1, 0, 'GrandGrandMean');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1,'gui','off'); % Make a new set
[ALLEEG,EEG, CURRENTSET]  = eeg_store(ALLEEG, EEG, CURRENTSET);                   % store it
Grand_Grand_Mean_Index    = CURRENTSET;                                           % keep track of it

eeglab redraw

%% Sorting things out over means and subjects
% First, the sequence of the two group means has be adjusted based on the grand grand mean
ALLEEG = pop_SortMSTemplates(ALLEEG, [Grand_Mean_AS_G_Index Grand_Mean_TA_G_Index], 1, Grand_Grand_Mean_Index);

% Then, we sort the individuals based on their group means
ALLEEG = pop_SortMSTemplates(ALLEEG, AS_Group_Index, 0, Grand_Mean_AS_G_Index); % AS Group
ALLEEG = pop_SortMSTemplates(ALLEEG, TA_Group_Index, 0, Grand_Mean_TA_G_Index); % TA Group

eeglab redraw

%% Savings
for i = 1:numel(ALLEEG)
    EEG = eeg_retrieve(ALLEEG,i);
    pop_saveset(EEG, 'filename', EEG.setname, 'filepath', SavePath, 'savemode', 'onefile');
end

%% Computing Parameters
for i = ClustPars.MinClasses:ClustPars.MaxClasses 
    FitPars.nClasses = i;
    % Using the the AS mean templates
    FileNameAS = ['AS_Mean_Template_Maps' num2str(i) '_Lambda' num2str(FitPars.lambda) '_b' num2str(FitPars.b) '.csv'];
    pop_QuantMSTemplates(ALLEEG, AS_Group_Index, 1, FitPars, numel(AllSubjects)+1, fullfile(SavePath,FileNameAS));

    % Using the the TA mean templates
    FileNameTA = ['TA_Mean_Template_Maps' num2str(i) '_Lambda' num2str(FitPars.lambda) '_b' num2str(FitPars.b) '.csv'];
    pop_QuantMSTemplates(ALLEEG, TA_Group_Index, 1, FitPars, numel(AllSubjects)+2, fullfile(SavePath,FileNameTA));

    % Using the grand grand mean template
    FileNameGGM = ['GG_Mean_Template_Maps' num2str(i) '_Lambda' num2str(FitPars.lambda) '_b' num2str(FitPars.b) '.csv'];
    pop_QuantMSTemplates(ALLEEG, AllSubjects, 1, FitPars, numel(AllSubjects)+3, fullfile(SavePath,FileNameGGM));
end