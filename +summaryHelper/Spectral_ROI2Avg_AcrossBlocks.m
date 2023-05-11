function [roi2avg] = Spectral_ROI2Avg_AcrossBlocks(dataTable, band)
%% Documentation
%
%   Create a ROI to average specified frequency band spectral power map
%       For each ROI, get all spectral power values in the specified 
%       frequency band in a segment. Then, average across that segment. 
%       Then, average across all segments in a block. 
%       Finally, store this value in the map 
%
%   ***THEN, average across all subjects with multiple blocks of the same
%       condition***
%
%   Make sure the data table includes the mapped condition names already 
%       --> make a table column called "mapCondition"
%

%% Matching Band to Frequencies
bandFreq = zeros(1,2);
if strcmp(band, 'Delta') == 1
    bandFreq = [0 4];
elseif strcmp(band, 'Theta') == 1
    bandFreq = [4 8];
elseif strcmp(band, 'Alpha') == 1
    bandFreq = [8 14];
elseif strcmp(band, 'Beta') == 1
    bandFreq = [14 30];
elseif strcmp(band, 'Gamma') == 1
    bandFreq = [30 50];
elseif strcmp(band, 'High Gamma') == 1
    bandFreq = [70 110];
elseif strcmp(band, 'All Gamma') == 1
    bandFreq = [30 120];
elseif strcmp(band, 'Total') == 1
    bandFreq = [0 120];
end

%% Creating the Contact2Average
% FOR EVERY SINGLE BLOCK IN THE DATA TABLE:
% 1st column is patient
% 2nd column is condition
% 3rd column is ROI to Average Band Spectral Power Map

[patientList, ia, ib] = unique(dataTable.patientID,'stable');

blockCells = cell(length(dataTable.patientID), 3);
blockCounter = 0;
for iPatient = 1:length(patientList)
    patientConditions = dataTable.mapCondition(ib==iPatient);
    for iCondition = 1:length(patientConditions)
        blockCounter = blockCounter + 1;
        currPatient = patientList{iPatient};
        currCondition = dataTable.mapCondition{blockCounter};
        currECoGchannels = dataTable.ECoGchannels{blockCounter};
        currSpecCell = dataTable.specAnalysis{blockCounter};
        specNames = fieldnames(currSpecCell);

        % Getting list of all ROIs
        ROIs = cell(1, length(currECoGchannels));
        for i = 1:length(currECoGchannels)
            ROIs{i} = currECoGchannels(i).oldROI;
        end

        % Creating the ROI to average band spectral power map
        %   Extract all values out of a segment, average across that
        %   segment, then average across all segments
        ROI2Avg = containers.Map();
        allROIs = unique(ROIs);
        for iROI = 1:length(allROIs)
            currROI = allROIs{iROI};
            rowIndices = [];
            for nROI = 1:length(ROIs)
                if strcmp(ROIs(nROI), currROI) == 1
                    rowIndices = [rowIndices nROI];
                end
            end
            
            allMeanValues = [];
            for iSegStruct = 1:size(specNames)
                currSegStruct = currSpecCell.(specNames{iSegStruct});
                            
                columnIndices = [];
                for i = 1:length(currSegStruct.freq)
                    if (currSegStruct.freq(i) >= bandFreq(1)) && (currSegStruct.freq(i) <= bandFreq(2))
                        columnIndices = [columnIndices i];
                    end
                end

                specValues = currSegStruct.powspctrm(rowIndices, columnIndices);
                allMeanValues = [allMeanValues mean(specValues(:), 'omitnan')];
            end
            ROI2Avg(currROI) = mean(allMeanValues(:), 'omitnan');
        %end of channels loop
        end

        blockCells{blockCounter, 1} = currPatient;
        blockCells{blockCounter, 2} = currCondition;
        blockCells{blockCounter, 3} = ROI2Avg;
    %end of conditions loop
    end
%end of patients loop
end

%% Create map of subject to conditions to be used below
% And, count how many unique combinations of subject/conditions there are

uniqueCounter = 0;
patient2conditions = containers.Map();
for iCell = 1:length(blockCells)
    % If the current patient is not in the map keys, add them and their
    % current condition
    currPatient = blockCells{iCell, 1};
    currCondition = string(blockCells{iCell, 2});
    if any(strcmp(currPatient, keys(patient2conditions))) == 0
        patient2conditions(currPatient) = currCondition;
        uniqueCounter = uniqueCounter + 1;
    else
        % If the current patient is in the map keys, check that whether the current condition is in there 
        % If not, add it to the array of conditions
        if any(strcmp(currCondition, patient2conditions(currPatient))) == 0
            currArray = [patient2conditions(currPatient) currCondition];
            patient2conditions(currPatient) = currArray;
            uniqueCounter = uniqueCounter + 1;
        end
    end
end

%% Combining subject blocks with same condition
% 1st column is subject
% 2nd column is condition
% 3rd column is ROI to avg across blocks map

comboCells = cell(uniqueCounter, 3);
comboCounter = 0;
patientList = unique(dataTable.patientID,'stable');
for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    currConditions = patient2conditions(currPatient);
    for iCondition = 1:length(currConditions)
        currCondition = currConditions(iCondition);

        % To store the array of all averages
        arrayROIMap = containers.Map();
        % To store the averages of the averages across subjects
        outputROIMap = containers.Map();

        % Find all cells with the current patient & condition combination
        for iCell = 1:length(blockCells)
            if strcmp(currPatient, blockCells{iCell, 1}) && strcmp(currCondition, blockCells{iCell, 2}) 
                currMap = blockCells{iCell, 3};
                allROIs = keys(currMap);
                for iROI = 1:length(allROIs)
                    currROI = allROIs{iROI};
                    currAvg = currMap(currROI);
                    
                    % If the ROI is already in arrayContactMap, add
                    % this average value to the array stored there
                    %
                    % Else, create the start of a new average array for 
                    % this ROI
                    if any(strcmp(currROI, keys(arrayROIMap))) == 1
                        currAvgArray = arrayROIMap(currROI);
                        currAvgArray = [currAvgArray currAvg];
                        arrayROIMap(currROI) = currAvgArray;
                    else
                        arrayROIMap(currROI) = currAvg;
                    end
                end
            end
        end
        % End of blockCells loop
        
        % Average across all block averages per ROI
        allROIs = keys(arrayROIMap);
        for iROI = 1:length(allROIs)
            currROI = allROIs{iROI};
            currAvgArray = arrayROIMap(currROI);
            outputROIMap(currROI) = mean(currAvgArray, 'omitnan');
        end
        
        comboCounter = comboCounter + 1;
        comboCells{comboCounter, 1} = currPatient;
        comboCells{comboCounter, 2} = currCondition;
        comboCells{comboCounter, 3} = outputROIMap;
    end
    % End of conditions loop
end
% End of patients loop


%% Converting it to a table to output
roi2avg = cell2table(comboCells, 'VariableNames', {'patientID', 'mapCondition', 'ROI2Avg'});

end