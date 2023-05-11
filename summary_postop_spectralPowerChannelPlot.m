function [dataTable] = summary_postop_spectralPowerChannelPlot(measure,DataSet,varargin)
%% Documentation
%
% Plots the average functional connectivity for a particular measure and 
% frequency across both conditions
%
% Averages each segment, across segments, across blocks for a 
% given subject & condition, and then across all channels
%
% varargin is the classic postHocOptionsParser parameters until the last
% pair which is 'Band'
%
%% Parsing varargin and loading data table

for i=1:length(varargin)
    if isequal('Band', varargin{i})
        band = varargin{i+1};
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

o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{1:end-4});
if ~iscell(o.DataSet), o.DataSet = {o.DataSet}; end

[dataTable] = loadGoodData(o);

[patientList,ia,ib] = unique(dataTable.patientID,'stable');

%% Mapping Condition Names

%Find index of 'StateMap' and mapping the original conditions onto the new ones
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

%% Creating the contact to spectral power average maps
% See functions below for more information

%contact2avg = summaryHelper.Spectral_Contact2Avg_PerBlock(dataTable, band);
contact2avg = summaryHelper.Spectral_Contact2Avg_AcrossBlocks(dataTable, band);

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

%% Plotting Conditions against each other

h = figure('Position',[65 158 1485 792]);

firstPlotBoolean = 0;
plotOrder = [];
allXPlotted = [];
allYPlotted = [];
for iPatient = 1:length(patientList)
    lateBoolean = 0;
    earlyBoolean = 0;
    dlrmBoolean = 0;
    currPatient = patientList{iPatient};
    for iCell = 1:length(contact2avg.patientID)
        if strcmp(currPatient, contact2avg.patientID{iCell}) == 1
            if contains(contact2avg.mapCondition{iCell}, 'Late') == 1
                lateBoolean = 1;
                firstDict = contact2avg.Contact2Avg{iCell};
                Channel2ROI = contact2avg.Contact2ROI{iCell};
            end
            if contains(contact2avg.mapCondition{iCell}, 'Early') == 1
                earlyBoolean = 1;
                secondDict = contact2avg.Contact2Avg{iCell};
                Channel2ROI = contact2avg.Contact2ROI{iCell};
            end
            if contains(contact2avg.mapCondition{iCell}, 'Dlrm') == 1
                dlrmBoolean = 1;
            end
        end
    end
    
    if earlyBoolean == 1 && lateBoolean == 1 
        if firstPlotBoolean == 0
            firstPlotBoolean = 1;
            if dlrmBoolean == 1
                plotOrder = ["Dlrm+", "Dlrm-"];
            elseif dlrmBoolean == 0
                plotOrder = ["Dlrm-", "Dlrm+"];
            end
        end

        currChannelsCell = keys(firstDict);

        figure(h)
        hold on;
        
        allLatePlotted = [];
        allEarlyPlotted = [];
        for iChannel = 1:length(currChannelsCell)
            currChannel = currChannelsCell{iChannel};
            if any( strcmp(currChannel, keys(secondDict)) ) == 0
                continue
            end

            currROI = Channel2ROI(currChannel);
            if strcmp(currROI, "SZ") == 1 || strcmp(currROI, "WM") == 1
                continue;
            end
            if any(strcmp(currROI, keys(colorMap))) == 0
                continue
            end
    
            firstSpectral = firstDict(currChannel);
            allLatePlotted = [allLatePlotted firstSpectral];
            
            secondSpectral = secondDict(currChannel);
            allEarlyPlotted = [allEarlyPlotted secondSpectral];
        end
        
        Late_average = mean(allLatePlotted, 'omitnan');
        Early_average = mean(allEarlyPlotted, 'omitnan');

        allXPlotted = [allXPlotted Late_average];
        allYPlotted = [allYPlotted Early_average];

        if dlrmBoolean == 1
            plot(Late_average, Early_average, '.', 'Color', [1 0 0], 'MarkerSize', 25);
        elseif dlrmBoolean == 0
            plot(Late_average, Early_average, '.', 'Color', [0 0 1], 'MarkerSize', 25);
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

xlabel("Late Average Spectral");
ylabel("Early Average Spectral");

legend(plotOrder);

plot([currMin currMax], [currMin currMax], 'g--', 'LineWidth', 1, 'DisplayName', 'y = x line');

if optionSetBoolean == 0
    sgtitle(['Average Subject ' band ' Band Spectral Power']);
elseif optionSetBoolean == 1
    sgtitle([optionSet ' Average Subject ' band ' Band Spectral Power']);
end


%% Saving Figure
if ~exist(['Z:\Ephys\AnalysisOutput\SpectralPower\Summaries\Channel Color\' DataSet], 'dir')
    mkdir(['Z:\Ephys\AnalysisOutput\SpectralPower\Summaries\Channel Color\' DataSet])
    error(['New Dataset Directory Created: ' DataSet ' --> Check if accidental before rerunning'])
end
fileOutPath = ['Z:\Ephys\AnalysisOutput\SpectralPower\Summaries\Channel Color\' DataSet '\summary_' band '_' optionSet '_spectralChannelPlot.fig'];
savefig(h, fileOutPath);
close(h);

end