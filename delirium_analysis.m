%% Getting IDs of all patients

dlrmMetaData = readtable("Z:\Ephys\MetaData\dataset\delirium.xlsx");
dlrmPatients = unique(dlrmMetaData.patientID);

ictalMetaData =  readtable("Z:\Ephys\MetaData\dataset\postIctal.xlsx");
ictalPatients = unique(ictalMetaData.patientID);

sleepMetaData = readtable("Z:\Ephys\MetaData\dataset\sleep.xlsx");
sleepPatients = unique(sleepMetaData.patientID);


combinedPatients = dlrmPatients;
for i = 1:length(ictalPatients)
    combinedPatients{end+1} = ictalPatients{i};
end
for i = 1:length(sleepPatients)
    combinedPatients{end+1} = sleepPatients{i};
end
combinedPatients = unique(combinedPatients);

%% Loading in subjects and blocks from two delirium datasets

postictal = readtable("Z:\Ephys\MetaData\dataset\postIctal.xlsx");
delirium = readtable("Z:\Ephys\MetaData\dataset\delirium.xlsx");

ictal_subjects = unique(postictal.patientID);
dlrm_subjects = unique(delirium.patientID);

all_subjects = ictal_subjects;
for i = 1:length(dlrm_subjects)
    if any(strcmp(dlrm_subjects(i), ictal_subjects)) == 1
        % no nothing
    else
        all_subjects{end+1} = dlrm_subjects{i};
    end
end

ictal_blocks = unique(postictal.block);
dlrm_blocks = unique(delirium.block);

all_blocks = ictal_blocks;
for i = 1:length(dlrm_blocks)
    if any(strcmp(dlrm_blocks(i), ictal_blocks)) == 1
        % no nothing
    else
        all_blocks{end+1} = dlrm_blocks{i};
    end
end

%% Initial Analysis Run
% runFuncs = {@wPLI_dbt_delta,@wPLI_dbt_theta,@wPLI_dbt_alpha,@wPLI_dbt_beta,...
%     @envCorrDBTorth_theta,@envCorrDBTorth_alpha,...
%     @envCorrDBTorth_beta,@envCorrDBTorth_gamma,@envCorrDBTorth_highGamma};
% runFuncs = {@wPLI_dbt_alpha,@envCorrDBTorth_gamma,@envCorrDBTorth_highGamma};

%runFuncs = {@specAnalysis,@wPLI_dbt_alpha,@envCorrDBTorth_gamma,@envCorrDBTorth_highGamma};
runFuncs = {@specAnalysis, @envCorrDBTorth_gamma};

runFuncs = {@specAnalysis, @wPLI_dbt_alpha, @wPLI_dbt_delta, @envCorrDBTorth_gamma, @envCorrDBTorth_highGamma};

runFuncs = {@wPLI_dbt_beta, @wPLI_dbt_theta, @wPLI_dbt_gamma, @envCorrDBTorth_theta, @envCorrDBTorth_alpha, @envCorrDBTorth_beta};
runFuncs = {@wPLI_dbt_delta, @envCorrDBTorth_theta, @envCorrDBTorth_alpha, @envCorrDBTorth_beta};


%patientAnalysis.runAnalysis(runFuncs,'Conditions',{'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate'},'OptionSet','segLen10');
%patientAnalysis.runAnalysis(runFuncs,'Conditions',{'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate''},'OptionSet','earlySVDsegLen10');
patientAnalysis.runAnalysis(runFuncs,'Conditions',{'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate'});
patientAnalysis.runAnalysis(runFuncs,'Conditions',{'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate'},'OptionSet','earlySVD');

%patientAnalysis.runAnalysis(runFuncs,'Conditions',{'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},'OptionSet','segLen10');
%patientAnalysis.runAnalysis(runFuncs,'Conditions',{'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'}},'OptionSet','earlySVDsegLen10');
patientAnalysis.runAnalysis(runFuncs,'Conditions',{'RS','RSresolved','RSctrlLate','RSictal','RSictalDlrm'});
patientAnalysis.runAnalysis(runFuncs,'Conditions',{'RS','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},'OptionSet','earlySVD');
% Took out RSpreEX because it's in the delirium one above

patientAnalysis.runAnalysis(runFuncs,'Subjects',combinedPatients);
patientAnalysis.runAnalysis(runFuncs,'Subjects',combinedPatients,'OptionSet','earlySVD');
%patientAnalysis.runAnalysis(runFuncs,'Subjects',combinedPatients,'OptionSet','preSVD');

% To test the runFuncs on just one subject:
patientAnalysis.runAnalysis(runFuncs,'Subjects', {'423L'}, 'Blocks', {'065'});
patientAnalysis.runAnalysis(runFuncs,'Subjects', {'423L'}, 'Blocks', {'065'}, 'OptionSet','earlySVD');

runFuncs = {@specAnalysis, @wPLI_dbt_alpha, @wPLI_dbt_delta, @envCorrDBTorth_gamma, @envCorrDBTorth_highGamma};
%runFuncs = {@wPLI_dbt_alpha, @envCorrDBTorth_gamma, @envCorrDBTorth_highGamma};

%patientAnalysis.runAnalysis(runFuncs,'Subjects', all_subjects, 'Blocks', all_blocks);
%patientAnalysis.runAnalysis(runFuncs,'Subjects', all_subjects, 'Blocks', all_blocks, 'OptionSet','earlySVD');

patientAnalysis.runAnalysis(runFuncs,'Subjects', {'439B'}, 'Blocks', all_blocks);
patientAnalysis.runAnalysis(runFuncs,'Subjects', {'439B'}, 'Blocks', all_blocks, 'OptionSet','earlySVD');

%% Post-op State Average
% Note that 439B's electrode reimplantation is dealth with in stateAverage
% for post-op

%yourStateMap should be a containers.Map so 
deliriumStateMap = containers.Map({'RSdlrm','RSresolved','RSctrlEarly','RSpreEx', 'RSctrlLate'},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
%run without any argument and it'll give you a file selector prompt
patientSummaries.stateAverage('envCorrDBTorth_highGamma','post-op delirium','StateMap',deliriumStateMap);
patientSummaries.stateAverage('envCorrDBTorth_gamma','post-op delirium','StateMap',deliriumStateMap);
patientSummaries.stateAverage('wPLI_dbt_alpha','post-op delirium','StateMap',deliriumStateMap);
patientSummaries.stateAverage('envCorrDBTorth_highGamma','post-op delirium','StateMap',deliriumStateMap,'OptionSet','earlySVD');
patientSummaries.stateAverage('envCorrDBTorth_gamma','post-op delirium','StateMap',deliriumStateMap,'OptionSet','earlySVD');
patientSummaries.stateAverage('wPLI_dbt_alpha','post-op delirium','StateMap',deliriumStateMap,'OptionSet','earlySVD');
%patientSummaries.stateAverage('envCorrDBTorth_highGamma','post-op delirium','StateMap',deliriumStateMap,'OptionSet','segLen10');
%patientSummaries.stateAverage('envCorrDBTorth_gamma','post-op delirium','StateMap',deliriumStateMap,'OptionSet','segLen10');
%patientSummaries.stateAverage('wPLI_dbt_alpha','post-op delirium','StateMap',deliriumStateMap,'OptionSet','segLen10');
%atientSummaries.stateAverage('envCorrDBTorth_highGamma','post-op delirium','StateMap',deliriumStateMap,'OptionSet','earlySVDsegLen10');
%patientSummaries.stateAverage('envCorrDBTorth_gamma','post-op delirium','StateMap',deliriumStateMap,'OptionSet','earlySVDsegLen10');
%patientSummaries.stateAverage('wPLI_dbt_alpha','post-op delirium','StateMap',deliriumStateMap,'OptionSet','earlySVDsegLen10');

%% Post-ictal State Average
% Note that 439B's electrode reimplantation is dealth with in stateAverage
% for post-op

%yourStateMap should be a containers.Map so 
% NOTE: we need to change this state map 10-22-2021
postictalStateMap = containers.Map({'RS','RSpreEx','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl','Ictal','IctalDlrm'});
%run without any argument and it'll give you a file selector prompt
patientSummaries.stateAverage('envCorrDBTorth_highGamma','post-ictal delirium','StateMap',postictalStateMap);
patientSummaries.stateAverage('envCorrDBTorth_gamma','post-ictal delirium','StateMap',postictalStateMap);
patientSummaries.stateAverage('wPLI_dbt_alpha','post-ictal delirium','StateMap',postictalStateMap);
patientSummaries.stateAverage('envCorrDBTorth_highGamma','post-ictal delirium','StateMap',postictalStateMap,'OptionSet','earlySVD');
patientSummaries.stateAverage('envCorrDBTorth_gamma','post-ictal delirium','StateMap',postictalStateMap,'OptionSet','earlySVD');
patientSummaries.stateAverage('wPLI_dbt_alpha','post-ictal delirium','StateMap',postictalStateMap,'OptionSet','earlySVD');
%patientSummaries.stateAverage('envCorrDBTorth_highGamma','post-ictal delirium','StateMap',postictalStateMap,'OptionSet','segLen10');
%patientSummaries.stateAverage('envCorrDBTorth_gamma','post-ictal delirium','StateMap',postictalStateMap,'OptionSet','segLen10');
%patientSummaries.stateAverage('wPLI_dbt_alpha','post-ictal delirium','StateMap',postictalStateMap,'OptionSet','segLen10');
%patientSummaries.stateAverage('envCorrDBTorth_highGamma','post-ictal delirium','StateMap',deliriumStateMap,'OptionSet','earlySVDsegLen10');
%patientSummaries.stateAverage('envCorrDBTorth_gamma','post-ictal delirium','StateMap',deliriumStateMap,'OptionSet','earlySVDsegLen10');
%patientSummaries.stateAverage('wPLI_dbt_alpha','post-ictal delirium','StateMap',deliriumStateMap,'OptionSet','earlySVDsegLen10');
%% Plotting Average Adjacency Matrix

% Naming for this function is all messed up with the new delirium dataset names
patientSummaries.plotAvgAdjMat; 

%
%
%% EFFECTIVE DIMENSIONALITY PLOTS
%
%
%
%
%
%
%% Effective Dimensionality Time Series within subjects across two conditions
% measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma'};
%
%measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha'};
%deliriumStateMap = containers.Map({'RSdlrm','RSresolved','RSctrlEarly','RSpreEx', 'RSctrlLate'},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'});
%for iMeas = 1:length(measureSet)
%    thisMeas = measureSet{iMeas};
%    patientSummaries.effDimTimeseries(thisMeas,'post-op delirium',[],'StateMap',deliriumStateMap,'OptionSet','earlySVDsegLen10');
%    patientSummaries.effDimTimeseries(thisMeas,'post-op delirium',[],'StateMap',deliriumStateMap,'OptionSet','segLen10');
%end
% 
% postictalStateMap = containers.Map({'RS','RSpreEx','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl','Ictal','IctalDlrm'});
% for iMeas = 1:length(measureSet)
%     thisMeas = measureSet{iMeas};
%     patientSummaries.effDimTimeseries(thisMeas,'post-ictal delirium',[],'StateMap',postictalStateMap);
% end
% 
% %% Effective Dimensionality Time Series for all with in one condition group
% % measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma'};
% measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha'};
% postictalStateMap = containers.Map({'RS','RSpreEx','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl','Ictal','IctalDlrm'});
% for iMeas = 1:length(measureSet)
%     thisMeas = measureSet{iMeas};
%     patientSummaries.effDimTimeseriesPerCondition(thisMeas,'post-ictal delirium',[],'StateMap',postictalStateMap);
% end
% 
% %% Mean Effective Dimensionality per Subject
% measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha'};
% postictalStateMap = containers.Map({'RS','RSpreEx','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl','Ictal','IctalDlrm'});
% for iMeas = 1:length(measureSet)
%     thisMeas = measureSet{iMeas};
%     patientSummaries.effDimPlotMean(thisMeas,'post-ictal delirium','StateMap',postictalStateMap);
% end
% 
% %% RS Effective Dimensionality vs Post-Ictal Effective Dimensionality (Color of Point = Delirious vs. Non-Delirious)
% measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha'};
% postictalStateMap = containers.Map({'RS','RSpreEx','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl','Ictal','IctalDlrm'});
% for iMeas = 1:length(measureSet)
%     thisMeas = measureSet{iMeas};
%     patientSummaries.effDimRSvsPostIctal(thisMeas,'post-ictal delirium','StateMap',postictalStateMap);
% end
%
%
%
%
%
%
%
%
%
% END EFFECTIVE DIMENSIONALITY PLOTS
%
%
%% SPECTRAL POWER and FUNCTIONAL CONNECTIVITY SUBPLOTS
%
%
%
%
%
%
%
%% Ictal vs IctalDlrm all within subject combinations spectral plot 
% 
% thisMeas = 'envCorrDBTorth_gamma';
% powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% 
% postictalStateMap = containers.Map({'RS','RSpreEx','RSctrlEarly','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
% 
% % First condition is x axis, second condition is y-axis
% for iBand = 1:length(powerBands)
%     thisBand = powerBands{iBand};
%     patientSummaries.spectralChannelPlot_ictaldlrm(thisMeas,'postIctal_combinations','StateMap', postictalStateMap, 'Band', thisBand);
% end
%
%% Generalizeable -- Plot ROI Spectral Power across chosen conditions

thisMeas = 'envCorrDBTorth_gamma';
powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

% ROI --> we want only intra-hemisphere spectral power!!!
% First condition is x axis, second condition is y-axis
for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Band', thisBand, 'Conditions', {'Ctrl', 'Ictal'});
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Band', thisBand, 'Conditions', {'Ctrl', 'IctalDlrm'});
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Band', thisBand, 'Conditions', {'Ictal', 'IctalDlrm'});
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Band', thisBand, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Band', thisBand, 'Conditions', {'DlrmLate', 'DlrmEarly'});

    patientSummaries.spectralPowerROIPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'Ctrl', 'Ictal'});
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'Ctrl', 'IctalDlrm'});
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'Ictal', 'IctalDlrm'});
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.spectralPowerROIPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'DlrmLate', 'DlrmEarly'});
end

%% Generalizeable -- Plot Channel Spectral Power across chosen conditions

thisMeas = 'envCorrDBTorth_gamma';
powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

% 'Hemisphere', 'both'

% First condition is x axis, second condition is y-axis
for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand, 'Conditions', {'Ctrl', 'Ictal'});
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand, 'Conditions', {'Ctrl', 'IctalDlrm'});
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand, 'Conditions', {'Ictal', 'IctalDlrm'});
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand, 'Conditions', {'DlrmLate', 'DlrmEarly'});

    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'Ctrl', 'Ictal'});
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'Ctrl', 'IctalDlrm'});
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'Ictal', 'IctalDlrm'});
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.spectralPowerChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'DlrmLate', 'DlrmEarly'});
end

%% Generalizeable -- Plot Spectral Power Ratio of Two Bands across chosen conditions

thisMeas = 'envCorrDBTorth_gamma';
% Choose any ratio of two of the following bands: {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% Put the slower band first in the order

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

% First condition is x axis, second condition is y-axis
    
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'Ctrl', 'Ictal'});
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'Ctrl', 'IctalDlrm'});
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'Ictal', 'IctalDlrm'});
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'CtrlLate', 'CtrlEarly'});
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'DlrmLate', 'DlrmEarly'});

patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'Ctrl', 'Ictal'});
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'Ctrl', 'IctalDlrm'});
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'Ictal', 'IctalDlrm'});
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'CtrlLate', 'CtrlEarly'});
patientSummaries.spectralBandRatio_ChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Bands', {'Delta', 'Beta'}, 'Conditions', {'DlrmLate', 'DlrmEarly'});

%% Generalizeable -- Plot Average Functional Connectivity across chosen conditions
% Can also do the average functional connectivity
% Generalizeable code

% OPTION FOR SVD DENOISING!!!
% Don't do both hemispheres for functional connectivity!!

%measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha', 'wPLI_dbt_delta'};
measureSet = {'wPLI_dbt_beta', 'wPLI_dbt_theta', 'wPLI_dbt_gamma', 'envCorrDBTorth_theta', 'envCorrDBTorth_alpha', 'envCorrDBTorth_beta'};
postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate'},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

% First condition is x axis, second condition is y-axis
for iMeas = 1:length(measureSet)
    thisMeas = measureSet{iMeas};
    patientSummaries.funcConnChannelPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'Conditions', {'Ctrl', 'Ictal'});
    patientSummaries.funcConnChannelPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'Conditions', {'Ctrl', 'IctalDlrm'});
    patientSummaries.funcConnChannelPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'Conditions', {'Ictal', 'IctalDlrm'});
    patientSummaries.funcConnChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.funcConnChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Conditions', {'DlrmLate', 'DlrmEarly'}); 

    patientSummaries.funcConnChannelPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'OptionSet','earlySVD', 'Conditions', {'Ctrl', 'Ictal'});
    patientSummaries.funcConnChannelPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'OptionSet','earlySVD', 'Conditions', {'Ctrl', 'IctalDlrm'});
    patientSummaries.funcConnChannelPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'OptionSet','earlySVD', 'Conditions', {'Ictal', 'IctalDlrm'});
    patientSummaries.funcConnChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'OptionSet','earlySVD', 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.funcConnChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'OptionSet','earlySVD', 'Conditions', {'DlrmLate', 'DlrmEarly'});    
end

%% Generalizeable -- Plot Average Functional Connectivity across chosen conditions
% Can also do the average functional connectivity
% Generalizeable code

% OPTION FOR SVD DENOISING!!!
% Don't do both hemispheres for functional connectivity!!

%measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha', 'wPLI_dbt_delta'};
measureSet = {'wPLI_dbt_beta', 'wPLI_dbt_theta', 'wPLI_dbt_gamma', 'envCorrDBTorth_theta', 'envCorrDBTorth_alpha', 'envCorrDBTorth_beta'};

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

% First condition is x axis, second condition is y-axis
for iMeas = 1:length(measureSet)
    thisMeas = measureSet{iMeas};
    patientSummaries.funcConnROIPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'Conditions', {'Ctrl', 'Ictal'});
    patientSummaries.funcConnROIPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'Conditions', {'Ctrl', 'IctalDlrm'});
    patientSummaries.funcConnROIPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'Conditions', {'Ictal', 'IctalDlrm'});
    patientSummaries.funcConnROIPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.funcConnROIPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Conditions', {'DlrmLate', 'DlrmEarly'}); 
    
    patientSummaries.funcConnROIPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'OptionSet','earlySVD', 'Conditions', {'Ctrl', 'Ictal'});
    patientSummaries.funcConnROIPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'OptionSet','earlySVD', 'Conditions', {'Ctrl', 'IctalDlrm'});
    patientSummaries.funcConnROIPlot(thisMeas,'post-ictal delirium','StateMap',  postictalStateMap, 'OptionSet','earlySVD', 'Conditions', {'Ictal', 'IctalDlrm'});
    patientSummaries.funcConnROIPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'OptionSet','earlySVD', 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.funcConnROIPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'OptionSet','earlySVD', 'Conditions', {'DlrmLate', 'DlrmEarly'});  
end

%% Generalizeable -- Plot Functional Connectivity of all ROIs to chosen ONE ROI across chosen conditions
% All subjects with correct conditions must have ROI
% If not the case, code will print out ROIs that all subjects share
% 
% measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha'};
% postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
% deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
% 
% % Need Two Conditions
% % Need ROI
% for iMeas = 1:length(measureSet)
%     thisMeas = measureSet{iMeas};
%     patientSummaries.funcConnSelectROI_AllSubj(thisMeas,'post-ictal delirium','StateMap',postictalStateMap,'Conditions',{'Ictal', 'IctalDlrm'},'ROI','InsA');
% end
% 
% %% Generalizeable -- Same as above but...only plots subjects that work
% % Below will plot all subjects that have the conditions and ROI
% % Doesn't extract any extra info like above, just plots these subjects
% 
% measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha'};
% postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
% deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
% 
% for iMeas = 1:length(measureSet)
%     thisMeas = measureSet{iMeas};
%     patientSummaries.funcConnSelectROI_SomeSubj(thisMeas,'post-ictal delirium','StateMap', postictalStateMap,'Conditions',{'Ictal', 'IctalDlrm'},'ROI','PMC');
% end
%
%% Spectral Ratio Plots -- Post-op -- early:late vs N3:WS for delirium + and -
% three delirium+: 423L, 439B, and 460L
% five delirium-: 403L, 409L, 418R, 567R, 585
% 
% thisMeas = 'envCorrDBTorth_gamma';
% powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% 
% deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
% % Not currently working
% for iBand = 1:length(powerBands)
%     thisBand = powerBands{iBand};
%     patientSummaries.postop_spectralRatioPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand);
% end
%
%% Post-op Spectral Ratio Plot V2 -- One condition at a time
% three delirium+: 423L, 439B, and 460L
% five delirium-: 403L, 409L, 418R, 567R, 585

% V2 is for the combinations of Early vs Wake/N3

thisMeas = 'envCorrDBTorth_gamma';
deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
early_conditions = {'DlrmEarly', 'CtrlEarly'};
%sleep_conditions = {'N3', 'WS'};
sleep_conditions = {'N3'};

for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    for iDlrm = 1:length(early_conditions)
        curr_early = early_conditions{iDlrm};
        for iSleep = 1:length(sleep_conditions)
            curr_sleep = sleep_conditions{iSleep};

            patientSummaries.postop_spectralRatioPlotV2(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', ...
                'Band', thisBand, 'EarlyCond', curr_early, "SleepCond", curr_sleep);
            patientSummaries.postop_spectralRatioPlotV2(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', ...
                'OptionSet', 'earlySVD', 'Band', thisBand, 'EarlyCond', curr_early, "SleepCond", curr_sleep);
        end
    end
end


%% Spectral Band Ratio State Ratio Plots -- Post-op -- early:late vs N3:WS for delirium + and -
% three delirium+: 423L, 439B, and 460L
% five delirium-: 403L, 409L, 418R, 567R, 585
% 
% thisMeas = 'envCorrDBTorth_gamma';
% %powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% 
% deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
% 
% patientSummaries.postop_spectralBandRatioRatioPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Bands', {'Delta', 'Beta'});
%
%% Spectral Ratio Plots -- Post-ictal -- ictal:ctrl vs N3:WS for delirium + and -
% two delirium+: 403L, 439B
% five delirium-: 418R, 423L, 514L, 567R, 585L
% 
% thisMeas = 'envCorrDBTorth_gamma';
% powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% 
% postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
% 
% for iBand = 1:length(powerBands)
%     thisBand = powerBands{iBand};
%     patientSummaries.postictal_spectralRatioPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand);
% end
%
%% Post-ictal Spectral Ratio Plot V2 -- One condition at a time
% two delirium+: 403L, 439B
% five delirium-: 418R, 423L, 514L, 567R, 585L

% V2 is for the combinations of Ictal/IctalDlrm vs Wake/N3

thisMeas = 'envCorrDBTorth_gamma';
postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});

powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
dlrm_conditions = {'Ictal', 'IctalDlrm'};
%sleep_conditions = {'N3', 'WS'};
sleep_conditions = {'N3'};

% ********** IN ITS CURRENT STATE WILL NOT GENERALIZE IF A SUBJECT IN BOTH
% THE SLEEP AND DELIRIUM DATASETS EXISTS THAT HAS BOTH A CTRL TO ICTAL AND
% A CTRL TO ICTALDLRM COMPARISON IT CAN MAKE **********
% --> will only make one of the two comparisons

for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    for iDlrm = 1:length(dlrm_conditions)
        curr_dlrm = dlrm_conditions{iDlrm};
        for iSleep = 1:length(sleep_conditions)
            curr_sleep = sleep_conditions{iSleep};

            patientSummaries.postictal_spectralRatioPlotV2(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', ...
                'Band', thisBand, 'DlrmCond', curr_dlrm, "SleepCond", curr_sleep);
            patientSummaries.postictal_spectralRatioPlotV2(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', ...
                'OptionSet', 'earlySVD', 'Band', thisBand, 'DlrmCond', curr_dlrm, "SleepCond", curr_sleep);
        end
    end
end

%% Spectral Band Ratio State Ratio Plots -- Post-ictal -- ictal:ctrl vs N3:WS for delirium + and -
% two delirium+: 403L, 439B
% five delirium-: 418R, 423L, 514L, 567R, 585L
% 
% thisMeas = 'envCorrDBTorth_gamma';
% powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% 
% postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
% 
% patientSummaries.postictal_spectralBandRatioRatioPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Bands', {'Delta', 'Beta'});
%
%% Average Spectral Power Time Series
% 
% thisMeas = 'envCorrDBTorth_gamma';
% powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% 
% postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
% deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
% 
% % Won't work with the current post-op dataset!!
% % Outputs every block of ictal dataset
% % Legend messed up...
% 
% for iBand = 1:length(powerBands)
%     thisBand = powerBands{iBand};
%     %patientSummaries.spectralTimeSeries(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand);
%     patientSummaries.spectralTimeSeries(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand);
% end
%
%
%
%
%
%
%
% END SPECTRAL POWER and FUNCTIONAL CONNECTIVITY SUBPLOTS
%
%
%% POWER SPECTRAL DENSITY PLOTS
%
%
%
%
%
%
%
%
%
%% Power Spectral Density Plot V3 -- ALL CHANNELS per subject

% Plot PSD on the vertical axis and frequency on the horizontal axis, both on log scales) 
% For both the delirium data and the sleep data
thisMeas = 'envCorrDBTorth_gamma';

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

patientSummaries.powerSpectralDensityV3(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both');
patientSummaries.powerSpectralDensityV3(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both');
patientSummaries.powerSpectralDensityV3(thisMeas,'sleep', 'Hemisphere', 'both');

patientSummaries.powerSpectralDensityV3(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD');
patientSummaries.powerSpectralDensityV3(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD');
patientSummaries.powerSpectralDensityV3(thisMeas,'sleep', 'Hemisphere', 'both', 'OptionSet', 'preSVD');
%
%
%
%
%
%
% END POWER SPECTRAL DENSITY PLOTS
%
%
%% SPECTRAL POWER AND FUNCTIONAL CONNECTIVITY SUMMARY PLOTS
%
%
%
%
%
%
%
%
%
%% Spectral Power BOX PLOTS V1
% Make sure no figures are currently open when running this!
% thisMeas = 'envCorrDBTorth_gamma';
% powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% 
% postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
% deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
% 
% % ***********V1 has each point be an entry in SpecAnalysis***********
% % UNSURE IF FUNCTIONAL RIGHT NOW
% for iBand = 1:length(powerBands)
%     thisBand = powerBands{iBand};
%     patientSummaries.spectralPower_BoxPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand);
%     patientSummaries.spectralPower_BoxPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand);
% end
%
%% Spectral Power BOX PLOTS V2
% Make sure no figures are currently open when running this!
thisMeas = 'envCorrDBTorth_gamma';
powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

% ***********V2 has each point be a channel average***********

for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    patientSummaries.spectralPower_BoxPlotV2(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand);
    patientSummaries.spectralPower_BoxPlotV2(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand);
    patientSummaries.spectralPower_BoxPlotV2(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand);
    patientSummaries.spectralPower_BoxPlotV2(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand);
end

%% Post-op Spectral Channel Plot Summary

thisMeas = 'envCorrDBTorth_gamma';
powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};

deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

% First condition is x axis, second condition is y-axis
for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    patientSummaries.summary_postop_spectralPowerChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.summary_postop_spectralPowerChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'CtrlLate', 'CtrlEarly'});
end

%% Post-ictal Spectral Channel Plot Summary

thisMeas = 'envCorrDBTorth_gamma';
powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});

% First condition is x axis, second condition is y-axis
for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    patientSummaries.summary_postictal_spectralPowerChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.summary_postictal_spectralPowerChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand, 'Conditions', {'CtrlLate', 'CtrlEarly'});
end

%% Post-op Functional Connectivity Plot Summary

measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha','wPLI_dbt_delta', 'wPLI_dbt_beta', 'wPLI_dbt_theta', 'wPLI_dbt_gamma', 'envCorrDBTorth_theta', 'envCorrDBTorth_alpha', 'envCorrDBTorth_beta'};

deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate',},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 

% First condition is x axis, second condition is y-axis
for iMeas = 1:length(measureSet)
    thisMeas = measureSet{iMeas};
    patientSummaries.summary_postop_funcConnChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.summary_postop_funcConnChannelPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'OptionSet','earlySVD', 'Conditions', {'CtrlLate', 'CtrlEarly'});
end


%% Post-ictal Functional Connectivity Plot Summary

measureSet = {'envCorrDBTorth_gamma','envCorrDBTorth_highGamma','wPLI_dbt_alpha', 'wPLI_dbt_delta'};
%measureSet = {'wPLI_dbt_beta', 'wPLI_dbt_theta', 'wPLI_dbt_gamma', 'envCorrDBTorth_theta', 'envCorrDBTorth_alpha', 'envCorrDBTorth_beta'};

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});

% First condition is x axis, second condition is y-axis
for iMeas = 1:length(measureSet)
    thisMeas = measureSet{iMeas};
    patientSummaries.summary_postictal_funcConnChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Conditions', {'CtrlLate', 'CtrlEarly'});
    patientSummaries.summary_postictal_funcConnChannelPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'OptionSet','earlySVD', 'Conditions', {'CtrlLate', 'CtrlEarly'});
end


%% Post-op Spectral Ratio Plot Summary
% 
% thisMeas = 'envCorrDBTorth_gamma';
% %powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% powerBands = {'Delta'};
% 
% deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate'},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
% 
% for iBand = 1:length(powerBands)
%     thisBand = powerBands{iBand};
%     patientSummaries.summary_postop_spectralRatioPlot(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand);
% end
% 
%
%% Post-ictal Spectral Ratio Plot Summary
% 
% thisMeas = 'envCorrDBTorth_gamma';
% %powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};
% powerBands = {'Delta'};
% 
% postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});
% 
% for iBand = 1:length(powerBands)
%     thisBand = powerBands{iBand};
%     patientSummaries.summary_postictal_spectralRatioPlot(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand);
% end
%
%% Post-op Spectral Ratio Plot Summary V2 -- ONLY Dlrm+/- vs N3!!!

thisMeas = 'envCorrDBTorth_gamma';

powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};

% Axis min and max are hard coded -- only work for Delta/Theta/Total

deliriumStateMap = containers.Map({'RSdlrm','RSresolved', 'RSctrlEarly','RSpreEx', 'RSctrlLate'},{'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate','CtrlLate'}); 
% Make sure to put Hemisphere before OptionSet in ordering
for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    patientSummaries.summary_postop_spectralRatioPlotV2(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'Band', thisBand)
    patientSummaries.summary_postop_spectralRatioPlotV2(thisMeas,'post-op delirium','StateMap', deliriumStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand)
end


%% Post-ictal Spectral Ratio Plot Summary V2 -- ONLY Dlrm+/- vs N3!!!

thisMeas = 'envCorrDBTorth_gamma';

powerBands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'High Gamma', 'All Gamma', 'Total'};

% Axis min and max are hard coded -- only work for Delta/Theta/Total

postictalStateMap = containers.Map({'RS','RSpreEx','RSresolved','RSctrlLate','RSictal','RSictalDlrm'},{'Ctrl','Ctrl','Ctrl', 'Ctrl', 'Ictal','IctalDlrm'});

for iBand = 1:length(powerBands)
    thisBand = powerBands{iBand};
    patientSummaries.summary_postictal_spectralRatioPlotV2(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'Band', thisBand)
    patientSummaries.summary_postictal_spectralRatioPlotV2(thisMeas,'post-ictal delirium','StateMap', postictalStateMap, 'Hemisphere', 'both', 'OptionSet', 'earlySVD', 'Band', thisBand)
end
%% Miscellaneous Debugging
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%% FIXING 439B ISSUE WITH REIMPLANTED ELECTRODES
% Checking whether the original and the reimplanted electrodes are within
% 5 mm of eachother in 3D space

% The CRU 005-053 Rev20221201BH -- Note BH not KN and the new date -- table has now been
% updated with this information
original_electrode_table = readtable("Z:\ElectrodeTables\439B_Electrode_Sites_KN.xlsx", "Sheet", "CRU 001-004 Rev20220801KN");
reimplanted_electrode_table = readtable("Z:\ElectrodeTables\439B_Electrode_Sites_KN.xlsx", "Sheet", "CRU 005-053 Rev20220801KN");

nrows = size(reimplanted_electrode_table);
nrows = nrows(1);

keep_array = [];
leave_array = [];
for i_row = 1:nrows
    if reimplanted_electrode_table.Contact(i_row) > 1000 
        if reimplanted_electrode_table.Channel(i_row) == original_electrode_table.Channel(i_row)

            original_x = original_electrode_table.RASX(i_row);
            original_y = original_electrode_table.RASY(i_row);
            original_z = original_electrode_table.RASZ(i_row);
    
            reimplanted_x = reimplanted_electrode_table.RASX(i_row);
            reimplanted_y = reimplanted_electrode_table.RASY(i_row);
            reimplanted_z = reimplanted_electrode_table.RASZ(i_row);
    
            euclidean_distance = sqrt( (reimplanted_x - original_x)^2 + (reimplanted_y - original_y)^2 + (reimplanted_z - original_z)^2 );

            disp([string(reimplanted_electrode_table.Contact(i_row)) ' Euclidean Distance: ' string(euclidean_distance)])

            if euclidean_distance < 5
                keep_array = [keep_array; reimplanted_electrode_table.Contact(i_row)];
            else
                leave_array = [leave_array; reimplanted_electrode_table.Contact(i_row)];
            end
        else
            disp([i_row ' -- row channels do not match'])
        end
    end
end

%% Now, of those that are in leave_array, find the closest channels to the new position


for n = 1:length(leave_array)
    curr_contact = leave_array(n);
    for i_row = 1:nrows 
        if reimplanted_electrode_table.Contact(i_row) == curr_contact
            min_euclidean_distance = 1000000000;

            reimplanted_x = reimplanted_electrode_table.RASX(i_row);
            reimplanted_y = reimplanted_electrode_table.RASY(i_row);
            reimplanted_z = reimplanted_electrode_table.RASZ(i_row);

            for j_row = 1:nrows
                original_x = original_electrode_table.RASX(j_row);
                original_y = original_electrode_table.RASY(j_row);
                original_z = original_electrode_table.RASZ(j_row);

                euclidean_distance = sqrt( (reimplanted_x - original_x)^2 + (reimplanted_y - original_y)^2 + (reimplanted_z - original_z)^2 );
                if euclidean_distance < min_euclidean_distance
                    min_euclidean_distance = euclidean_distance;
                    closest_contact = original_electrode_table.Contact(j_row);
                end
            end
        end
    end

    disp(append(string(curr_contact), "'s closest contact in first block is ", string(closest_contact), " with a Euclidean Distance of ", string(min_euclidean_distance)))
end

%% Checking which sides are dominant now
blockone_electrode_table = readtable("Z:\ElectrodeTables\439B_Electrode_Sites_KN.xlsx", "Sheet", "CRU 001-004 Rev20220801KN");
blocktwo_electrode_table = readtable("Z:\ElectrodeTables\439B_Electrode_Sites_KN.xlsx", "Sheet", "CRU 005-053 Rev20221201BH");

nrows = size(blocktwo_electrode_table);
nrows = nrows(1);

blockone_R_counter = 0;
blockone_L_counter = 0;
blocktwo_R_counter = 0;
blocktwo_L_counter = 0;

for i_row = 1:nrows
    % If the contact number is less than 1000, include it in the counts
    if blocktwo_electrode_table.Contact(i_row) < 1000 
        if strcmp(blockone_electrode_table.Side(i_row), "R") == 1
            blockone_R_counter = blockone_R_counter + 1;
        elseif strcmp(blockone_electrode_table.Side(i_row), "L") == 1
            blockone_L_counter = blockone_L_counter + 1;
        end

        if strcmp(blocktwo_electrode_table.Side(i_row), "R") == 1
            blocktwo_R_counter = blocktwo_R_counter + 1;
        elseif strcmp(blocktwo_electrode_table.Side(i_row), "L") == 1
            blocktwo_L_counter = blocktwo_L_counter + 1;
        end
    end
end
disp(' ')
disp(['Block one has ' char(string(blockone_R_counter)) ' Right Side Channels'])
disp(['Block one has ' char(string(blockone_L_counter)) ' Left Side Channels'])
disp(['Block two has ' char(string(blocktwo_R_counter)) ' Right Side Channels'])
disp(['Block two has ' char(string(blocktwo_L_counter)) ' Left Side Channels'])

