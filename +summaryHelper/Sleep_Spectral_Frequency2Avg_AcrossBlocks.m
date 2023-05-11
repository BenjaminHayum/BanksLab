function [freq2avg] = Sleep_Spectral_Frequency2Avg_AcrossBlocks(dataTable)
%% Documentation
%
%   Create a frequency to average spectral power map
%       For all ROIs that aren't SZ or WM and are in the ROI Color Keys, 
%       get all spectral values for the frequency in that segment.
%       Then, average across that segment. 
%       Then, average across all segments in a block. 
%       Finally, store this value in the map 
%
%   ***THEN, average across all subjects with multiple blocks of the same
%       condition*** 
%
%   Make sure the data table includes the mapped condition names already 
%       --> make a table column called "mapCondition"
%
%   
%   ** If you want to edit this to do a Power Spectral Density for only
%   select ROI, edit lines 60-65 to only get row indices from what you
%   choose **
%
%% Creating the ColorMap to use to filter out channel indices
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
%% Creating the Frequency2Average
% FOR EVERY SINGLE BLOCK IN THE DATA TABLE:
% 1st column is Patient
% 2nd column is Condition
% 3rd column is Frequency to Average Spectral Power

[patientList, ia, ib] = unique(dataTable.patientID,'stable');

blockCells = cell(length(dataTable.patientID)*5, 2);
blockCounter = 0;
for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    patientRows = ib==iPatient;

    currSpecCell = dataTable.specAnalysis{iPatient};
    specNames = fieldnames(currSpecCell);

    % Only including channels that aren't in SZ or WM ROIs and are ROIs
    % That are in the ROI Color Map (as has been standard)
    patientECoGchannels = dataTable.ECoGchannels{iPatient};
    goodRowIndices = [];
    for iROI = 1:length(patientECoGchannels)
        currROI = patientECoGchannels(iROI).oldROI;
        if strcmp(currROI, "SZ") == 0 && strcmp(currROI, "WM") == 0 && any(strcmp(currROI, keys(ROIcolorMap))) == 1
            goodRowIndices = [goodRowIndices iROI];
        end
    end

    % Listing the states for this patient
    patientStates = dataTable.states(patientRows);
    patientStates = patientStates{:};
    allStates = unique(patientStates);
    for iCondition = 1:length(allStates)
        currState = allStates{iCondition};
        blockCounter = blockCounter + 1;
        
        % Initializing the Maps (to be used in DictionaryCells)
        Freq2Avg = containers.Map();
        
        if strcmp(currState, 'WS') == 1
            stateKey = 'W_';
        elseif strcmp(currState, 'N1') == 1
            stateKey = 'N1_';
        elseif strcmp(currState, 'N2') == 1
            stateKey = 'N2_';
        elseif strcmp(currState, 'N3') == 1
            stateKey = 'N3_';
        elseif strcmp(currState, 'R') == 1
            stateKey = 'R_';
        end
        
        % Creating the frequency to average spectral power map
        %   Extract all values out of a segment, average across that
        %   segment, then average across all segments
        frequencies = currSpecCell.(specNames{1}).freq;
        for iFreq = 1:length(frequencies)
            currFreq = num2str(frequencies(iFreq));
            columnIndex = iFreq;

            allMeanSegValues = [];
            currStateSegCounter = -1;
            for iSegStruct = 1:size(specNames)
                currSegName = specNames{iSegStruct};
                if contains(currSegName, stateKey) == 1
                    currStateSegCounter = currStateSegCounter + 1;
                    currSegStruct = currSpecCell.(append(stateKey, num2str(currStateSegCounter)));
                    
                    % The column index is just  the frequency index
                    specValues = currSegStruct.powspctrm(goodRowIndices, columnIndex);
                    allMeanSegValues = [allMeanSegValues mean(specValues(:), 'omitnan')];
                end
            end
            % Making a map of each ROI and its corresponding mean spectral
            % power
            Freq2Avg(currFreq) = mean(allMeanSegValues(:), 'omitnan');
        end
           
        % Inputting it all
        blockCells{blockCounter, 1} = string(currPatient);
        blockCells{blockCounter, 2} = string(currState);            
        blockCells{blockCounter, 3} = Freq2Avg;
    %End of states loop
    end  
%End of patients loop     
end

%% Create map of subject to conditions to be used below
% And, count how many unique combinations of subject/conditions there are

uniqueCounter = 0;
patient2conditions = containers.Map(); % Goes unused in this code
for iCell = 1:length(blockCells)
    if isempty(blockCells{iCell, 1}) == 1
        break
    end
    % If the current patient is not in the map keys, add them and their
    % current condition
    currPatient = string(blockCells{iCell, 1});
    currState = string(blockCells{iCell, 2});
    if any(strcmp(currPatient, keys(patient2conditions))) == 0
        patient2conditions(currPatient) = currState;
        uniqueCounter = uniqueCounter + 1;
    else
        % If the current patient is in the map keys, check that whether the current condition is in there 
        % If not, add it to the array of conditions
        if any(strcmp(currState, patient2conditions(currPatient))) == 0
            currArray = [patient2conditions(currPatient) currState];
            patient2conditions(currPatient) = currArray;
            uniqueCounter = uniqueCounter + 1;
        end
    end
end

%% Combining subject blocks with same condition
% 1st column is Patient
% 2nd column is Condition
% 3rd column is Frequency to Average Spectral Power

comboCells = cell(uniqueCounter, 3);
% For each patient and condition:
comboCellCounter = 0;
for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    patientRows = ib==iPatient;
    patientStates = dataTable.states(patientRows);
    patientStates = patientStates{:};
    patientConditions = unique(patientStates);
    for iCondition = 1:length(patientConditions)
        currState = patientConditions{iCondition};

        % To store the array of all averages
        arrayFreqMap = containers.Map();
        % To store the averages of the averages across subjects
        outputFreqMap = containers.Map();

        % Find all cells with the current patient & condition combination
        for iCell = 1:length(blockCells)
            if strcmp(currPatient, blockCells{iCell, 1}) == 1 && strcmp(currState, blockCells{iCell, 2}) == 1
                currMap = blockCells{iCell, 3};
                allFreqs = keys(currMap);
                for iFreq = 1:length(allFreqs)
                    currFreq = num2str(allFreqs{iFreq});
                    currAvg = currMap(currFreq);

                    % If the frequency is already in arrayFreqMap, add
                    % this average value to the array stored there
                    %
                    % Else, create the start of a new average array for 
                    % this frequency
                    if any(strcmp(currFreq, keys(arrayFreqMap))) == 1
                        currAvgArray = arrayFreqMap(currFreq);
                        currAvgArray = [currAvgArray currAvg];
                        arrayFreqMap(currFreq) = currAvgArray;
                    else
                        arrayFreqMap(currFreq) = currAvg;
                    end
                end

            end
        end
        %End of dictionaryCells loop
        
        % Average across all block averages per frequency
        allFreqs = keys(arrayFreqMap);
        for iFreq = 1:length(allFreqs)
            currFreq = num2str(allFreqs{iFreq});
            currAvgArray = arrayFreqMap(currFreq);
            outputFreqMap(currFreq) = mean(currAvgArray, 'omitnan');
        end
        
        comboCellCounter = comboCellCounter + 1;
        comboCells{comboCellCounter, 1} = currPatient;
        comboCells{comboCellCounter, 2} = currState;
        comboCells{comboCellCounter, 3} = outputFreqMap;

    end
    %End of conditions loop
end
% End of patients loop

%% Output

freq2avg = cell2table(comboCells, 'VariableNames', {'patientID', 'mapCondition', 'Freq2Avg'});

end