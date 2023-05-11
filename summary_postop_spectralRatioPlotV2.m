function [dlrm_dataTable] = summary_postop_spectralRatioPlotV2(measure,DataSet,varargin)
%% Documentation
%
% Find subjects that are in both datatables
% Label dlrm+ vs dlrm-
% Get the average spectral power from early vs late and N3/N2 vs WS
%   Dictionary with all 4 values as well as the subject name and
%   delirium+ or delirium-
%
% ****Plot the Dlrm+/Dlrm- against N3****
%
% Plots the average spectral power for a particular band averaged across
% every contact/channel
%
% Averages each segment, across segments, and across blocks for a 
% given subject & condition.
%
%% Parsing varargin
for i=1:length(varargin)
    if isequal('Band', varargin{i})
        band = varargin{i + 1};
    end
end

% Checking if there's an OptionSet
optionSetBoolean = 0;
for i=1:length(varargin)
    if isequal('OptionSet', varargin{i})
        optionSet = varargin{i+1};
        optionSetBoolean = 1;
        break
    end
end
if optionSetBoolean == 0
    optionSet = 'NoOptionSet';
end
%% Loading Delirium Datatable
o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{1:end-2});

[dlrm_dataTable] = loadGoodData(o);

[dlrm_PatientList,dlrm_ia,dlrm_ib] = unique(dlrm_dataTable.patientID,'stable');

%% Loading Sleep Datatable
hemisphereBoolean = 0;
for iArg = 1:length(varargin)
    if strcmp(varargin{iArg}, 'Hemisphere') == 1
        hemisphereBoolean = 1;
        hemisphere = varargin{iArg + 1};
        break
    end
end

if hemisphereBoolean == 0 && optionSetBoolean == 0
    o = ecogutils.postHocOptionsParser(measure, 'sleep');
elseif hemisphereBoolean == 1 && optionSetBoolean == 0
    o = ecogutils.postHocOptionsParser(measure, 'sleep', 'Hemisphere', hemisphere);
elseif hemisphereBoolean == 0 && optionSetBoolean == 1
    o = ecogutils.postHocOptionsParser(measure, 'sleep', 'OptionSet', optionSet);
elseif hemisphereBoolean == 1 && optionSetBoolean == 1
    if strcmp(optionSet, 'earlySVD') == 1
        o = ecogutils.postHocOptionsParser(measure, 'sleep', 'Hemisphere', hemisphere, 'OptionSet', 'preSVD');
        optionSet = 'earlySVD-preSVD';
    else
        o = ecogutils.postHocOptionsParser(measure, 'sleep', 'Hemisphere', hemisphere, 'OptionSet', optionSet);
    end
end

[sleep_dataTable] = loadGoodData(o);

[sleep_PatientList, sleep_ia,sleep_ib] = unique(sleep_dataTable.patientID,'stable');

%% Finding the subjects that both datasets have
shared_PatientList = [];
for iPatient = 1:length(dlrm_PatientList)
    if any(strcmp(dlrm_PatientList{iPatient}, sleep_PatientList)) == 1
        shared_PatientList = [shared_PatientList string(dlrm_PatientList{iPatient})];
    end
end

%% Mapping Conditions
%Find index of 'StateMap' and mapping the original conditions onto the new ones
stateMapIndex = -1;
for i=1:length(varargin)
    if isequal('StateMap', varargin{i})
        stateMapIndex = i;
        break
    end
end

dlrm_dataTable.mapCondition = cell(length(dlrm_dataTable.conditions), 1);
if stateMapIndex ~= -1
    for i = 1:length(dlrm_dataTable.conditions)
        dlrm_dataTable.mapCondition{i} = varargin{stateMapIndex+1}(dlrm_dataTable.conditions{i});
    end
else
    for i = 1:length(dlrm_dataTable.conditions)
        dlrm_dataTable.mapCondition{i} = dlrm_dataTable.conditions{i};
    end
end

%% All the possible ROIs to iterate through

T = readtable('Z:\Ephys\MetaData\ROITable\ROI_table.xlsx');
oldROIArray = T.oldROI;
colorR = T.colorR;
colorG = T.colorG;
colorB = T.colorB;

colorMap = containers.Map();
for iROI = 1:length(oldROIArray)
    currROI = oldROIArray{iROI};
    currRed = colorR(iROI);
    currGreen = colorG(iROI);
    currBlue = colorB(iROI);

    colorMap(currROI) = [currRed, currGreen, currBlue];
end

%% Loading delirium and sleep contact2avg's

dlrm_contact2avg = summaryHelper.Spectral_Contact2Avg_AcrossBlocks(dlrm_dataTable, band);
sleep_contact2avg = summaryHelper.Sleep_Spectral_Contact2Avg_AcrossBlocks(sleep_dataTable, band);

%% Putting it into the format that the plot code wants

% Column 1 is subject number
% Column 2 is dlrm+ vs dlrm-
% Column 3 is early band average spectral power Channel dictionary
% Column 4 is late band average spectral power Channel dictionary
% Column 5 is the channel to ROI dictionary
DLRMCells = cell(length(shared_PatientList), 5);

% Column 1 is Subject
% Column 2 is whether subject is N2 or N3
% Column 3 is N3/N2 band average spectral power Channel dictionary
% Column 4 is WS band average spectral power Channel dictionary
SLEEPCells = cell(length(shared_PatientList), 4);

% Removing patients that aren't in both datasets
for iPatient = 1:length(shared_PatientList)
    currPatient = shared_PatientList{iPatient};
    % Delirium Half
    for iDlrm = 1:length(dlrm_contact2avg.patientID)
        if strcmp(currPatient, dlrm_contact2avg.patientID{iDlrm}) == 1
            DLRMCells{iPatient, 1} = currPatient;
            if contains(dlrm_contact2avg.mapCondition{iDlrm}, "Ctrl") == 1
                DLRMCells{iPatient, 2} = "Dlrm -";
            elseif contains(dlrm_contact2avg.mapCondition{iDlrm}, "Dlrm") == 1
                DLRMCells{iPatient, 2} = "Dlrm +";
            end
            if contains(dlrm_contact2avg.mapCondition{iDlrm}, "Early") == 1
                DLRMCells{iPatient, 3} = dlrm_contact2avg.Contact2Avg{iDlrm};
            elseif contains(dlrm_contact2avg.mapCondition{iDlrm}, "Late") == 1
                DLRMCells{iPatient, 4} = dlrm_contact2avg.Contact2Avg{iDlrm};
            end
            DLRMCells{iPatient, 5} = dlrm_contact2avg.Contact2ROI{iDlrm};
        end
    end

    % Sleep Half
    for iSleep = 1:length(sleep_contact2avg.patientID)
        if strcmp(currPatient, sleep_contact2avg.patientID{iSleep}) == 1
            SLEEPCells{iPatient, 1} = currPatient;
            SLEEPCells{iPatient, 2} = sleep_contact2avg.state{iSleep};
            SLEEPCells{iPatient, 3} = sleep_contact2avg.N3N2contact2avg{iSleep};
            SLEEPCells{iPatient, 4} = sleep_contact2avg.WScontact2avg{iSleep};
        end
    end
end


%% Plotting Ratios against each other (color = dlrm+ vs dlrm-)

h = figure('Position',[65 158 1485 792]);

firstPlotBoolean = 0;
plotOrder = [];

firstNoDlrmDone = 0;
firstDlrmDone = 0;

allXPlotted = [];
allYPlotted = [];

% DlrmCells -- plot average of all matching channels dlrm+ and dlrm- color
% coded to N2/N3
for iPatient = 1:length(shared_PatientList)
    figure(h)
    hold on;

    currPatient = DLRMCells{iPatient, 1};
    
    % Delirium Status for determining the color of the plot
    dlrmStatus = DLRMCells{iPatient, 2};
    if strcmp(dlrmStatus, "Dlrm -") == 1
        dlrmBoolean = 0;
        if firstPlotBoolean == 0
            firstPlotBoolean = 1;
            plotOrder = {"Dlrm-", "Dlrm+"};
        end
    elseif strcmp(dlrmStatus, "Dlrm +") == 1
        dlrmBoolean = 1;
        if firstPlotBoolean == 0
            firstPlotBoolean = 1;
            plotOrder = {"Dlrm+", "Dlrm-"};
        end
    end
    
    dlrmEarlyMap = DLRMCells{iPatient, 3};
    sleepLateMap = SLEEPCells{iPatient, 3};

    dlrmEarlyKeys = keys(dlrmEarlyMap);
    sleepLateKeys = keys(sleepLateMap);
    
    Channel2ROI = DLRMCells{iPatient, 5};

    % Go through each ROI one by one and plot
    % Only plot the ROIs that both have
    allSubjSleep = [];
    allSubjDelirium = [];
    for iChannel = 1:length(dlrmEarlyKeys)
        currChannel = dlrmEarlyKeys{iChannel};
        % Checking if all dictionaries have the Channel
        plotBoolean = 0;
        if any(strcmp(currChannel, sleepLateKeys)) == 1
            plotBoolean = 1;

            currROI = Channel2ROI(currChannel);
            if any(strcmp(currROI, "SZ")) == 1 || strcmp(currROI, "WM") == 1
                continue
            end
            if any(strcmp(currROI, keys(colorMap))) == 0
                continue
            end
        end
        
        if plotBoolean == 1
            earlyAvg = dlrmEarlyMap(currChannel);
            N3Avg = sleepLateMap(currChannel);
                    
            % To later plot average of
            allSubjSleep = [allSubjSleep N3Avg];
            allSubjDelirium = [allSubjDelirium earlyAvg];
        end
    end
    %End of Channels Loop
    
    X_average = mean(allSubjSleep, 'omitnan');
    Y_average = mean(allSubjDelirium, 'omitnan');
    
    allXPlotted = [allXPlotted X_average];
    allYPlotted = [allYPlotted Y_average];
    % Hard-coded Min and Max
    %currMin = 10^4;
    %currMax = 10^6;

    
    if dlrmBoolean == 1
        if firstDlrmDone == 0
            plot(X_average, Y_average, '.', 'Color', [1 0 0], 'MarkerSize', 25, 'DisplayName', 'Dlrm+');
            firstDlrmDone = 1;
        else
            plot(X_average, Y_average, '.', 'Color', [1 0 0], 'MarkerSize', 25, 'HandleVisibility', 'off');
        end
    elseif dlrmBoolean == 0
        if firstNoDlrmDone == 0
            plot(X_average, Y_average, '.', 'Color', [0 0 1],'MarkerSize', 25, 'DisplayName', 'Dlrm-');
            firstNoDlrmDone = 1;
        else
            plot(X_average, Y_average, '.', 'Color', [0 0 1],'MarkerSize', 25, 'HandleVisibility', 'off');
        end
    end
end

[currMin, currMax] = summaryHelper.squareLogAxisLimits(allXPlotted,allYPlotted);

currMin = currMin/2;
currMax = currMax*2;

xlim([currMin currMax]);
set(gca, 'XScale', 'log');
ylim([currMin currMax]);
set(gca, 'YScale', 'log');
axis square;
%set(gca, 'XTick', [10^4, 10^5, 10^6])
%set(gca, 'XMinorTick', 'on')
%set(gca, 'YTick', [10^4, 10^5, 10^6])
%set(gca, 'YMinorTick', 'on')

ylabel("Early Average Spectral Power");
xlabel("N3/N2 Average Spectral Power");

legend('show');

if optionSetBoolean == 0
    sgtitle(['postop ' band ' Band Subject Average Channel Spectral Power']);
elseif optionSetBoolean == 1
    sgtitle(['postop ' optionSet ' ' band ' Band Subject Average Channel Spectral Power']);
end


plot([currMin currMax], [currMin currMax], 'g--', 'LineWidth', 1, 'DisplayName', 'y = x line');

%% Saving Figure

if ~exist(['Z:\Ephys\AnalysisOutput\SpectralPower\Summaries\Sleep Comparison\Channel Color\' DataSet], 'dir')
    mkdir(['Z:\Ephys\AnalysisOutput\SpectralPower\Summaries\Sleep Comparison\Channel Color\' DataSet])
    error(['New Dataset Directory Created: ' DataSet ' --> Check if accidental before rerunning'])
end
fileOutPath = ['Z:\Ephys\AnalysisOutput\SpectralPower\Summaries\Sleep Comparison\Channel Color\' DataSet '\summary_' band '_' optionSet '_sleepChannelComparison.fig'];
savefig(h, fileOutPath);
close(h);


end