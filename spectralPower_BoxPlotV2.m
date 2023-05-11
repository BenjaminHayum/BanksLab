function [dataTable] = spectralPower_BoxPlotV2(measure,DataSet,varargin)
%% Documentation
%
% Create a box plot where each point is a channel average
% Order boxplots by median channel average spectral power
%
% Do Wilcoxin Rank Sum test at the end on ordering
%

%% Parsing varargin and loading data table
for i=1:length(varargin)
    if isequal('Band', varargin{i})
        band = varargin{i+1};
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

o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{1:end-2});
if ~iscell(o.DataSet), o.DataSet = {o.DataSet}; end

[dataTable] = loadGoodData(o);

[patientList,ia,ib] = unique(dataTable.patientID,'stable');

%% Mapping Condition Names
% Find index of 'StateMap'
stateMapIndex = -1;
for i=1:length(varargin)
    if isequal('StateMap', varargin{i})
        stateMapIndex = i;
        break
    end
end
if strcmp(DataSet, 'sleep') == 0
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

conditionList = keys(colorMap);

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

%% Creating the contact to spectral power average maps
% See functions below for more information

contact2avg = summaryHelper.Spectral_Contact2Avg_PerBlock(dataTable, band);
%contact2avg = summaryHelper.Spectral_Contact2Avg_AcrossBlocks(dataTable, band);

%% Combining subject blocks in same condition

% 1st column is subject
% 2nd column is condition
% 3rd column is ROI avg across blocks

prePlotCells = cell(length(patientList)*length(conditionList), 3);
% For each patient and condition:
plotCellCounter = 0;
for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    for iCondition = 1:length(conditionList)
        currCondition = conditionList{iCondition};

        % To store the array of all averages
        arrayChannelMap = containers.Map();
        % To store the averages of the averages across subjects
        outputChannelMap = containers.Map();
        Channel2ROI_V2 = containers.Map();

        % Go through all the cells and for the matching subject/condition
        % combos, put all their ROIs
        for iCell = 1:length(contact2avg.patientID)
            if strcmp(currPatient, contact2avg.patientID{iCell}) == 1 && strcmp(currCondition, contact2avg.mapCondition{iCell}) == 1
                
                % Go through all ROIs which should be the same within a
                % subject
                currMap = contact2avg.Contact2Avg{iCell};
                Channel2ROI = contact2avg.Contact2ROI{iCell};
                allChannels = keys(currMap);
                for iChannel = 1:length(allChannels)
                    currChannel = allChannels{iChannel};
                    currROI = Channel2ROI(currChannel);
                    if strcmp(currROI, "SZ") == 1 || strcmp(currROI, "WM") == 1
                        continue;
                    end
                    if any(strcmp(currROI, keys(ROIcolorMap))) == 0
                        continue
                    end
                    Channel2ROI_V2(currChannel) = currROI;

                    currAvg = currMap(currChannel);
                    if any(strcmp(currChannel, keys(arrayChannelMap))) == 1
                        currAvgArray = arrayChannelMap(currChannel);
                        currAvgArray = [currAvgArray currAvg];
                        arrayChannelMap(currChannel) = currAvgArray;
                    else
                        arrayChannelMap(currChannel) = currAvg;
                    end
                end
            end
        end
        %End of contact2avg loop
        
        % Calculating averages and assigning to plotCells
        if isempty(keys(arrayChannelMap)) == 0
            allChannels = keys(arrayChannelMap);
            for iChannel = 1:length(allChannels)
                currChannel = allChannels{iChannel};
                if any(strcmp(currChannel, keys(Channel2ROI_V2))) == 1
                    currAvgArray = arrayChannelMap(currChannel);
                    outputChannelMap(currChannel) = mean(currAvgArray, 'omitnan');
                end
            end
        end

        plotCellCounter = plotCellCounter + 1;
        prePlotCells{plotCellCounter, 1} = currPatient;
        prePlotCells{plotCellCounter, 2} = currCondition;
        %prePlotCells{plotCounter, 3} = values(outputROIMap); Need array below
        cellValues = values(outputChannelMap);
        arrayValues = zeros(1, length(cellValues));
        for iVal = 1:length(cellValues)
            arrayValues(iVal) = cellValues{iVal};
        end
        prePlotCells{plotCellCounter, 3} = arrayValues;

    end
    %End of conditions loop
end
% End of patients loop

%% Creating new patient list sorted by highest median

% Sort via highest control median spectral power
if strcmp(DataSet, "post-op delirium") == 1
    ctrlCondition = "Late";
elseif strcmp(DataSet, "post-ictal delirium") == 1
    ctrlCondition = "Ctrl";
end

medianPatientList = cell(length(patientList), 1);

medianMap = containers.Map();
for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    for iCell = 1:length(prePlotCells)
        if strcmp(currPatient, prePlotCells{iCell, 1}) == 1 
            if isempty(prePlotCells{iCell, 3}) == 0 && contains(prePlotCells{iCell, 2}, ctrlCondition) == 1
                currMedian = median(prePlotCells{iCell, 3}, 'all', 'omitnan');
            end
        end
    end
    medianMap(currPatient) = currMedian;
end

sortedMedianVals = sort(cell2mat(values(medianMap)), 'descend');
for iVal = 1:length(sortedMedianVals)
    currMedian = sortedMedianVals(iVal);
    for iPatient = 1:length(patientList)
        currPatient = patientList{iPatient};
        if currMedian == medianMap(currPatient) 
            medianPatientList{iVal} = currPatient;
        end
    end
end


%% Turning it directly into what's to be plotted -- cell array of matrices
% Matrix set ups:
%   1 matrix per condition and each condition has 
%       length(patientList) columns 
%       Number of data points rows

plotCells = cell(1, length(conditionList));
% For axes limits
minVal = Inf;
maxVal = 0;

% Finding what the dimensions of each array are to be
matrixColumnNum = length(medianPatientList);
matrixRowNumMap = containers.Map();
for iCell = 1:length(contact2avg.patientID)
    currCondition = contact2avg.mapCondition{iCell};
    currLength = length(contact2avg.Contact2Avg{iCell});
    if ~any(strcmp(keys(matrixRowNumMap), currCondition)) == 1
        matrixRowNumMap(currCondition) = currLength;
    else
        if currLength > matrixRowNumMap(currCondition)
            matrixRowNumMap(currCondition) = currLength;
        end
    end
end

%Pre-allocating space for matrices
for iCell = 1:length(plotCells)
    plotCells{iCell} = NaN(matrixRowNumMap(conditionList{iCell}), matrixColumnNum);
end

%Inputting the values
for iCell = 1:length(prePlotCells)
    %Getting patient number
    currPatient = prePlotCells{iCell, 1};
    for iPatient = 1:length(medianPatientList)
        if strcmp(medianPatientList{iPatient}, currPatient) == 1
            patientNum = iPatient;
        end
    end
    %Getting condition number
    currCondition = prePlotCells{iCell, 2};
    for iCondition = 1:length(conditionList)
        if strcmp(conditionList{iCondition}, currCondition) == 1
            conditionNum = iCondition;
        end
    end

    %inputting it into plotCells
    currValues = prePlotCells{iCell, 3};
    for iVal = 1:length(currValues)
        % For axes limits
        if log10(currValues(iVal)) < minVal
            minVal = log10(currValues(iVal));
        end
        if log10(currValues(iVal)) > maxVal
            maxVal = log10(currValues(iVal));
        end
        plotCells{conditionNum}(iVal, patientNum) = log10(currValues(iVal));
    end
end


%% Plotting 

% boxplotGroup(x) receives a 1xm cell array where each element is a matrix with
% n columns and produces n groups of boxplot boxes with m boxes per group.

% Matrix set ups:
%   1 matrix per condition
%   Each condition has length(patientList) columns and number of data
%   points rows

% Secondary labels is patients
% Primary labels is conditions
% interGroupSpace to separate it

% Wants colors in a specific format:
colors = ones(length(conditionList), 3);
for iCondition = 1:length(conditionList)
    currTrio = colorMap(conditionList{iCondition});
    for iTrio = 1:length(currTrio)
        colors(iCondition, iTrio) = currTrio(iTrio);
    end
end

% Plotting!
h = boxplotGroup(plotCells, 'primaryLabels', conditionList, 'secondaryLabels', medianPatientList, 'groupLine', true, ...
    'Colors', colors, 'Width', 0.75, 'GroupType','betweenGroups', 'BoxStyle', 'filled', 'OutlierSize', 4, 'Symbol', '.');

h.axis.XTickLabelRotation = 90;
ylim(h.axis2, [minVal - 1, maxVal + 1])
% Adjusting the x-axis's many titles
try
    gf = gcf;
    gf.Children(1).XRuler.TickLabelGapOffset = 1;
    gf.Children(2).XRuler.TickLabelGapOffset = 35;
    gf.Children(2).XLabel.String = 'Subjects';
    gf.Children(1).XLabel.String = '';
catch
    error("Close all other figures before running :) ")
end

ylabel("Log(Channel Average Spectral Power)");
if optionSetBoolean == 0
    sgtitle([DataSet ' ' band ' Band Channel Average Spectral Power Boxplot']);
elseif optionSetBoolean == 1
    sgtitle([DataSet ' ' optionSet ' ' band ' Band Channel Average Spectral Power Boxplot']);
end


%% Saving the Plots

if ~exist(['Z:\Ephys\AnalysisOutput\SpectralPower\Boxplot\' DataSet], 'dir')
    mkdir(['Z:\Ephys\AnalysisOutput\SpectralPower\Boxplot\' DataSet])
    error(['New Dataset Directory Created: ' DataSet ' --> Check if accidental before rerunning'])
end
fileOutPath = ['Z:\Ephys\AnalysisOutput\SpectralPower\Boxplot\' DataSet '\' band '_' optionSet '_SpectralBoxplotV2.fig'];
savefig(h.figure, fileOutPath);
close(h.figure);

%% Doing a Wilcoxin Rank Sum Test on the Ranks of the Control Median Log Spectral Powers for the groups of Dlrm+ vs Dlrm-

% Sort via highest control median log(spectral power)
if strcmp(DataSet, "post-op delirium") == 1
    ctrlCondition = "Late";
elseif strcmp(DataSet, "post-ictal delirium") == 1
    ctrlCondition = "Ctrl";
end

NonDlrmVals = [];
DlrmVals = [];

for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    for iCell = 1:length(prePlotCells)
        if strcmp(currPatient, prePlotCells{iCell, 1}) == 1 
            if isempty(prePlotCells{iCell, 3}) == 0 && contains(prePlotCells{iCell, 2}, ctrlCondition) == 1
                currMedian = median(prePlotCells{iCell, 3}, 'all', 'omitnan');


                for jCell = 1:length(prePlotCells)
                    if strcmp(currPatient, prePlotCells{jCell, 1}) == 1 && contains(prePlotCells{jCell, 2}, "Dlrm") == 1
                        if isempty(prePlotCells{jCell, 3}) == 0
                            DlrmVals = [DlrmVals currMedian];
                        else
                            NonDlrmVals = [NonDlrmVals currMedian];
                        end
                    end
                end


            end
        end
    end
end

disp(' ')
disp("p-value for a two-sided Wilcoxon Rank Sum Test of the Control Log(Spectral Powers) for the groups of Dlrm+ vs Dlrm-:")
disp(' ')
disp(append("Dataset: ", DataSet))
disp(append("Band: ", band))
p_value = ranksum(DlrmVals, NonDlrmVals)

x = "Just for setting Break Points :)";

end