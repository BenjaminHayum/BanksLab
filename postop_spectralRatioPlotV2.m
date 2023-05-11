function [dlrm_dataTable] = postop_spectralRatioPlotV2(measure,DataSet,varargin)
%% Documentation
% Find subjects that are in both datatables
% Label dlrm+ vs dlrm-
% Get the average spectral power from early vs late and N3/N2 vs WS
    % Dictionary with all 4 values as well as the subject name and
    % delirium+ or delirium-
% Plot the specified conditions against each other
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
        break
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

for i=1:length(varargin)
    if isequal('EarlyCond', varargin{i})
        earlyCondition = varargin{i+1};
        break
    end
end

for i=1:length(varargin)
    if isequal('SleepCond', varargin{i})
        sleepCondition = varargin{i+1};
        break
    end
end

%% Loading Delirium Datatable
o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{1:end-6});

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

%Getting the subplotsSize and indices needed
subplotsSize = 0;
indexArray = [];
if strcmp(earlyCondition, 'DlrmEarly') == 1
    dlrmCounter = 0;
    for iCell = 1:length(DLRMCells)
        if strcmp(DLRMCells{iCell, 2}, "Dlrm +")
            dlrmCounter = dlrmCounter + 1;
            indexArray = [indexArray iCell];
        end
    end
    for iSqrt = 1:dlrmCounter
        squareSize = iSqrt*iSqrt;
        if squareSize >= dlrmCounter
            subplotsSize = iSqrt;
            break;
        end
    end
elseif strcmp(earlyCondition, 'CtrlEarly') == 1
    ctrlCounter = 0;
    for iCell = 1:length(DLRMCells)
        if strcmp(DLRMCells{iCell, 2}, "Dlrm -")
            ctrlCounter = ctrlCounter + 1;
            indexArray = [indexArray iCell];
        end
    end
    for iSqrt = 1:ctrlCounter
        squareSize = iSqrt*iSqrt;
        if squareSize >= ctrlCounter
            subplotsSize = iSqrt;
            break;
        end
    end
end

finalDlrmCells = DLRMCells(indexArray, :);
finalSleepCells = SLEEPCells(indexArray, :);


% Use indices below
for iPatient = 1:length(indexArray)
    figure(h)
    subplot(subplotsSize,subplotsSize, iPatient);
    hold on;

    currPatient = finalDlrmCells{iPatient, 1};
    
    % Delirium Status for determining the color of the plot
    earlyMap = finalDlrmCells{iPatient, 3};

    if strcmp(sleepCondition, "N3") == 1
        sleepMap = finalSleepCells{iPatient, 3};
    elseif strcmp(sleepCondition, "WS") == 1
        sleepMap = finalSleepCells{iPatient, 4};
    end

    allXPlotted = [];
    allYPlotted = [];

    earlyKeys = keys(earlyMap);
    sleepKeys = keys(sleepMap);

    Channel2ROI = finalDlrmCells{iPatient, 5};

    % Go through each ROI one by one and plot
    % Only plot the ROIs that both have
    for iChannel = 1:length(earlyKeys)
        currChannel = earlyKeys{iChannel};
        % Checking if all dictionaries have the Channel
        plotBoolean = 0;
        if any(strcmp(currChannel, sleepKeys)) == 1
            plotBoolean = 1;
    
            currROI = Channel2ROI(currChannel);
            if strcmp(currROI, "SZ") == 1 || strcmp(currROI, "WM") == 1
                continue;
            end
            if any(strcmp(currROI, keys(colorMap))) == 1
                currColor = colorMap(currROI);
            else
                continue
            end
        end
        
        if plotBoolean == 1
            earlyAvg = earlyMap(currChannel);
            sleepAvg = sleepMap(currChannel);
                    
            plot(sleepAvg,  earlyAvg, '.', 'Color', currColor, 'MarkerSize', 10.5); 

            % To later create the axes limits
            allXPlotted = [allXPlotted sleepAvg];
            allYPlotted = [allYPlotted earlyAvg];
        end
    end
    %End of Channels Loop

    [currMin, currMax] = summaryHelper.squareLogAxisLimits(allXPlotted,allYPlotted);

    xlim([currMin currMax]);
    set(gca, 'XScale', 'log');
    ylim([currMin currMax]);
    set(gca, 'YScale', 'log');
    if strcmp(finalDlrmCells{iPatient, 2}, "Dlrm +") == 1
        text(currMax, currMax, 'Dlrm+', 'Color', [1 0 0]);
        plot([currMin currMax], [currMin currMax], 'r--', 'LineWidth', 1);
    elseif strcmp(finalDlrmCells{iPatient, 2}, "Dlrm -") == 1
        text(currMax, currMax, 'Dlrm-', 'Color', [0 0 1]);
        plot([currMin currMax], [currMin currMax], 'b--', 'LineWidth', 1);
    end

    axis square;
    %set(gca, 'XTick', [10^3, 10^4, 10^5, 10^6])
    %set(gca, 'XMinorTick', 'on')
    %set(gca, 'YTick', [10^3, 10^4, 10^5, 10^6])
    %set(gca, 'YMinorTick', 'on')
    
    ylabel(append(earlyCondition, " Spectral Power"))
    if strcmp(sleepCondition, "N3") == 1
        xlabel(append(finalSleepCells{iPatient, 2}, " Spectral Power"));
    elseif strcmp(sleepCondition, "WS") == 1
        xlabel("WS Spectral Power");
    end

    title(currPatient);
    if optionSetBoolean == 0
        sgtitle([band ' Band Average Channel Spectral Power Sleep vs Early Post-op Comparison']);
    elseif optionSetBoolean == 1
        sgtitle([optionSet ' ' band ' Band Average Channel Spectral Power Sleep vs Early Post-op Comparison']);
    end
    
end


%% Saving Figure
if ~exist(['Z:\Ephys\AnalysisOutput\SpectralPower\Sleep Comparison\Channel Color\' DataSet], 'dir')
    mkdir(['Z:\Ephys\AnalysisOutput\SpectralPower\Sleep Comparison\Channel Color\' DataSet])
    error(['New Dataset Directory Created: ' DataSet ' --> Check if accidental before rerunning'])
end
fileOutPath = ['Z:\Ephys\AnalysisOutput\SpectralPower\Sleep Comparison\Channel Color\' DataSet '\' band '_' optionSet '_' earlyCondition 'vs' sleepCondition 'SpectralChannel.fig'];
savefig(h, fileOutPath);
close(h);

end