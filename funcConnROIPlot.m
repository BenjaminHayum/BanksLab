function [dataTable] = funcConnROIPlot(measure,DataSet,varargin)
%% Documentation
%
% Plots the average functional connectivity for a particular measure and frequency
% across every ROI and both conditions
%
% Averages each segment, across segments, and across blocks for a 
% given subject & condition.
%
% varargin is the classic postHocOptionsParser parameters until the last
% pair which is 'Condition'
%
%% Load data table and ROI table

o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{1:end-2});
if ~iscell(o.DataSet), o.DataSet = {o.DataSet}; end

[dataTable] = loadGoodData(o);

[patientList,ia,ib] = unique(dataTable.patientID,'stable');

%% Mapping Conditions and Finding the FC Frequency

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

%Getting Conditions
for i=1:length(varargin)
    if isequal('Conditions', varargin{i})
        mainConditions = varargin{i + 1};
        break
    end
end

%Finding the FC Frequency
if contains(o.measure,'wPLI_')
    freqs = fieldnames(dataTable.(o.measure){ia(1)}{1}.wPLI_debias);
    if length(freqs)>1, error('Only 1 frequency supported.'), end
    freq = freqs{1,1};
elseif contains(o.measure,'envCorrDBT')
    freqs = fieldnames(dataTable.(o.measure){ia(1)}{1}.envCorr);
    if length(freqs)>1, error('Only 1 frequency supported.'), end
    freq = freqs{1,1};
else
    error('Unknown/unsupported connMeasure');
end

%% Hard coding in 439B with same hemisphere to the dataset for Post-op Runs
if strcmp(o.DataSet, 'post-op delirium') == 1
    varargin{end-1} = 'Hemisphere';
    varargin{end} = 'L';
    new_o = ecogutils.postHocOptionsParser(measure,'439B',varargin{:});
    [new_table] = loadGoodData(new_o);

    for iRow = 1:length(dataTable.patientID)
        if strcmp(dataTable.patientID{iRow}, new_table.patientID{1}) && strcmp(dataTable.block{iRow}, new_table.block{1}) 
            dataTable(iRow, :) = new_table(1, :);
        end
        if strcmp(dataTable.patientID{iRow}, new_table.patientID{2}) && strcmp(dataTable.block{iRow}, new_table.block{2}) 
            dataTable(iRow, :) = new_table(2, :);
        end
    end
end

%% Find index of 'StateMap' and mapping the original conditions onto the new ones
stateMapIndex = -1;
for i=1:length(varargin)
    if isequal('StateMap', varargin{i})
        stateMapIndex = i;

        break
    end
end

dataTable.mapCondition = cell(length(dataTable.conditions), 1);
if stateMapIndex ~= -1
    for i = 1:length(dataTable.conditions)
        dataTable.mapCondition{i} = varargin{stateMapIndex+1}(dataTable.conditions{i});
    end
else
    for i = 1:length(dataTable.conditions)
        dataTable.mapCondition{i} = dataTable.conditions{i};
    end
end

%% Creating the ROI to functional connectivity average maps
% See functions below for more information

%roi2avg = summaryHelper.FuncConn_ROI2Avg_PerBlock(dataTable, o.measure, freq);
roi2avg = summaryHelper.FuncConn_ROI2Avg_AcrossBlocks(dataTable, o.measure, freq);

if contains(o.measure,'wPLI_')
    measure = 'wPLI';
elseif contains(o.measure, 'envCorr')
    measure = 'envCorr';
end

%% Creating the ColorMap

%Use oldROI and colorR/colorG/colorB in MetaData/ROI_table.xlsx

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

%% Plotting Conditions against Each Other

h = figure('Position',[65 158 1485 792]);

% Finding how many patients are used
usedPatientCounter = 0;
% Loop through each patient and check whether they have the specified
% condiitons. If so, add to the counter to be used to determine the
% subplotsize
for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    hasFirstCondition = false;
    hasSecondCondition = false;
    for iCell = 1:length(roi2avg.mapCondition)
        if strcmp(currPatient, roi2avg.patientID{iCell}) && strcmp(mainConditions{1}, roi2avg.mapCondition{iCell})
            hasFirstCondition = true;
        end
        if strcmp(currPatient, roi2avg.patientID{iCell}) && strcmp(mainConditions{2}, roi2avg.mapCondition{iCell})
            hasSecondCondition = true;
        end
    end
    if hasFirstCondition && hasSecondCondition
        usedPatientCounter = usedPatientCounter + 1;
    end
end
% Getting the subplotsSize
subplotsSize = 0;
for iSqrt = 1:usedPatientCounter
    squareSize = iSqrt*iSqrt;
    if squareSize >= usedPatientCounter
        subplotsSize = iSqrt;
        break;
    end
end

plotCounter = 0;
for iPatient = 1:length(patientList)
    firstBoolean = 0;
    secondBoolean = 0;
    currPatient = patientList{iPatient};
    
    % For each patient with the two conditions, extract their roi2avg maps
    for iCell = 1:length(roi2avg.patientID)
        if strcmp(currPatient, roi2avg.patientID{iCell}) == 1
            if strcmp(mainConditions{1}, roi2avg.mapCondition{iCell}) == 1
                firstBoolean = 1;
                firstDict = roi2avg.ROI2Avg{iCell};
            end
            if strcmp(mainConditions{2}, roi2avg.mapCondition{iCell}) == 1
                secondBoolean = 1;
                secondDict = roi2avg.ROI2Avg{iCell};
            end
        
        end
    end
  
    
    if firstBoolean == 1 && secondBoolean == 1
        currROIsCell = keys(firstDict);
        %Have to go index by index to convert the cell of strings into an array
        currROIsArray = string(zeros(1, length(currROIsCell)));
        for i = 1:length(currROIsCell)
            currROIsArray(1, i) = string(currROIsCell{1, i});
        end

        plotCounter = plotCounter + 1;
        figure(h)
        subplot(subplotsSize,subplotsSize, plotCounter);
        hold on;
        
        allFirstPlotted = [];
        allSecondPlotted = [];
        for iROI = 1:length(currROIsArray)
            currROI = currROIsArray(iROI);

            % Plot all ROIs that are in both roi2avg dictionaries,
            % aren't the ROIs of SZ or WM, and are in the ROI Color Map
            if any( strcmp(currROI, keys(secondDict)) ) == 0
                continue
            end
            if strcmp(currROI, "SZ") == 1 || strcmp(currROI, "WM") == 1
                continue;
            end
            if any(strcmp(keys(colorMap), currROI))
                currColor = colorMap(currROI);
    
                firstSpectral = firstDict(currROI);
                allFirstPlotted = [allFirstPlotted firstSpectral];
                
                secondSpectral = secondDict(currROI);
                allSecondPlotted = [allSecondPlotted secondSpectral];
    
                text(firstSpectral, secondSpectral, currROI, 'Color', currColor, 'FontSize', 7, 'HorizontalAlignment','center');
            end
        end
      
        [currMin, currMax] = summaryHelper.squareLogAxisLimits(allFirstPlotted,allSecondPlotted);

        xlim([currMin currMax]);
        set(gca, 'XScale', 'log');
        ylim([currMin currMax]);
        set(gca, 'YScale', 'log');
    
        axis square;
        plot([currMin currMax], [currMin currMax], 'b--', 'LineWidth', 1);
        xlabel(append(mainConditions{1}, " Average FC"));
        ylabel(append(mainConditions{2}, " Average FC"));
        title(patientList(iPatient));
        
        if optionSetBoolean == 0
            sgtitle([measure freq]);
        elseif optionSetBoolean == 1
            sgtitle([measure freq ' ' optionSet]);
        end
    end
end

legend("y = x Line");

%% Saving Plots
if ~exist(['Z:\Ephys\AnalysisOutput\FunctionalConnectivity\ROI Color\' DataSet], 'dir')
    mkdir(['Z:\Ephys\AnalysisOutput\FunctionalConnectivity\ROI Color\' DataSet])
    error(['New Dataset Directory Created: ' DataSet ' --> Check if accidental before rerunning'])
end
fileOutPath = ['Z:\Ephys\AnalysisOutput\FunctionalConnectivity\ROI Color\' DataSet '\' measure freq '_' optionSet '_'  upper(mainConditions{1}) 'vs' upper(mainConditions{2})  '.fig'];
savefig(h, fileOutPath);
close(h);

end