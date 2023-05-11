function [dataTable] = powerSpectralDensityV3(measure,DataSet,varargin)
%% Documentation
%
% Plots the average spectral power for each band frequency value across every
% contact/channel for every condition
%
% Averages each segment, across segments, and across blocks for a 
% given subject & condition.
%
% If you want to only plot PSDs for a specific ROI, edit the code in 
% "Spectral_Frequency2Avg_AcrossBlocks" or
% "Sleep_Spectral_Frequency2Avg_AcrossBlocks"
%
%% Parsing varargin and loading data table
%Find index of 'StateMap'
stateMapIndex = -1;
for i=1:length(varargin)
    if isequal('StateMap', varargin{i})
        stateMapIndex = i;
        break
    end
end

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

if stateMapIndex == -1
    o = ecogutils.postHocOptionsParser(measure,DataSet);
else
    o = ecogutils.postHocOptionsParser(measure,DataSet, varargin{1:2});
end

[dataTable] = loadGoodData(o);

[patientList,ia,ib] = unique(dataTable.patientID,'stable');

%% Mapping Condition Names
if strcmp(DataSet, 'sleep') == 0
    dataTable.mapCondition = cell(length(dataTable.conditions), 1);
    if stateMapIndex ~= -1
        for i = 1:length(dataTable.conditions)
            dataTable.mapCondition{i} = varargin{stateMapIndex+1}(dataTable.conditions{i});
        end
        conditionList = unique(dataTable.mapCondition, 'stable');
    else
        for i = 1:length(dataTable.conditions)
            dataTable.mapCondition{i} = dataTable.conditions{i};
        end
        conditionList = unique(dataTable.mapCondition, 'stable');
    end
end

%% Getting Colors for Dataset
switch DataSet
    case 'sleep'
        colorMap = ...
            containers.Map({'WS','N1','N2','N3','R'},...
            {[139,0,0]./255,[0,128,0]./255,[0,0,255]./255,[128,128,128]./255,[255,165,0]./255});
    case 'post-op delirium'
        colorMap = ...
            containers.Map({'DlrmEarly','DlrmLate','CtrlEarly','CtrlLate'},...
            {[0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840] } );
    case 'post-ictal delirium'
        colorMap = ...
            containers.Map({'Ctrl', 'Ictal', 'IctalDlrm'},...
            {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250] } );
    otherwise
        error('Need to add this DataSet to state/color map list')
end

%% Creating the ColorMap
%Use oldROI and colorR/colorG/colorB in MetaData/ROI_table.xlsx

T = readtable('Z:\Ephys\MetaData\ROITable\ROI_table.xlsx');
oldROIArray = T.oldROI;
colorR = T.colorR;
colorG = T.colorG;
colorB = T.colorB;

ROIcolorMap = containers.Map();
for iROI = 1:length(oldROIArray)
    currROI = oldROIArray{iROI};
    currRed = colorR(iROI);
    currGreen = colorG(iROI);
    currBlue = colorB(iROI);

    ROIcolorMap(currROI) = [currRed, currGreen, currBlue];
end

%% Creating the frequency to spectral power average maps
% See functions below for more information

if strcmp(DataSet, "sleep") == 1
    freq2avg = summaryHelper.Sleep_Spectral_Frequency2Avg_AcrossBlocks(dataTable);
else
    freq2avg = summaryHelper.Spectral_Frequency2Avg_AcrossBlocks(dataTable);
end

%% Plotting Conditions against each other

h = figure('Position',[65 158 1485 792]);

%Getting the subplotsSize
subplotsSize = 0;
for iSqrt = 1:length(patientList)
    squareSize = iSqrt*iSqrt;
    if squareSize >= length(patientList)
        subplotsSize = iSqrt;
        break;
    end
end

% Go through each cell and line up the patient name extracted against
% patientList
plotCounter = 0;
endArray = []; % For making the end legend correct
for iCell = 1:length(freq2avg.patientID)
    currPatient = char(freq2avg.patientID{iCell});
    currCondition = char(freq2avg.mapCondition{iCell});
    
    for iPatient = 1:length(patientList)
        if strcmp(patientList{iPatient}, currPatient) == 1
            plotCounter = iPatient;
        end
    end
    if plotCounter == length(patientList)
        endArray = [endArray string(currCondition)];
    end

    figure(h)
    subplot(subplotsSize,subplotsSize, plotCounter);
    hold on;

    currMap = freq2avg.Freq2Avg{iCell};
    tempFreq = keys(currMap);

    frequencies = zeros(1,length(tempFreq));
    for iFreq = 1:length(tempFreq)
        frequencies(iFreq) = str2num(tempFreq{iFreq});
    end
    frequencies = sort(frequencies, 'ascend');

    % Add all frequencies and their average spectral power to arrays to
    % be later plotted
    freqArray = [];
    specArray = [];
    for iFreq = 1:length(frequencies)
        currFreq = frequencies(iFreq);
        currSpec = currMap(num2str(currFreq));

        freqArray = [freqArray currFreq];
        specArray = [specArray currSpec];
    end
    
    currColor = colorMap(currCondition);
    plot(freqArray, specArray, 'Color', currColor);
    %plot(freqArray, specArray, 'Color', currColor, 'DisplayName', currCondition);
    
    xlim([0 125]);
    set(gca, 'XScale', 'log');
    ylim([0.5 5000000]);
    set(gca, 'YScale', 'log');
    set(gca, 'YTick', [1, 100, 10000, 1000000])
    set(gca, 'YMinorTick', 'on')

    xlabel("Frequencies");
    ylabel("Average Spectral Power");
    title(currPatient);
    
    if optionSetBoolean == 0
        sgtitle([DataSet ' Average Spectral Power Density Plot']); 
    elseif optionSetBoolean == 1
        sgtitle([DataSet ' ' optionSet ' Average Spectral Power Density Plot']); 
    end
end

% Custom thing at the end to make the legend
allColors = keys(colorMap);
allVals = values(colorMap);
legendArray = endArray;
for iColor = 1:length(allColors)
    currCondition = allColors{iColor};
    if any(strcmp(currCondition, endArray)) == 1
        continue;
    else
        plot(0, 0, 'Color', allVals{iColor});
        legendArray = [legendArray currCondition];
    end
end
legend(legendArray);
%legend();

%% Saving the Plots
if ~exist(['Z:\Ephys\AnalysisOutput\PowerSpectralDensity\' DataSet], 'dir')
    mkdir(['Z:\Ephys\AnalysisOutput\PowerSpectralDensity\' DataSet])
    error(['New Dataset Directory Created: ' DataSet ' --> Check if accidental before rerunning'])
end
fileOutPath = ['Z:\Ephys\AnalysisOutput\PowerSpectralDensity\' DataSet '\' DataSet '_' optionSet '_PSD_v3.fig'];
savefig(h, fileOutPath);
close(h);

end