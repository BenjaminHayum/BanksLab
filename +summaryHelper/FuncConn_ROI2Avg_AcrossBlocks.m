function [roi2avg] = FuncConn_ROI2Avg_AcrossBlocks(dataTable, measure, freq)
%% Documentation
%
%   Create an ROI to average specified measure and frequency functional 
%       connectivity map
%
%   For each ROI, get all functional connectivity values in a segment. 
%       Then, average across that segment. 
%       Then, average across all segments in a block. 
%       Finally, store this value in the map 
%
%   Have accomodating contact to ROI map for later color coded plots   
%
%   Make sure the data table includes the mapped condition names already 
%       --> make a table column called "mapCondition"
%


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
        currConnCell = dataTable.(measure){blockCounter};

        % Getting list of all ROIs
        ROIs = cell(1, length(currECoGchannels));
        for i = 1:length(currECoGchannels)
            ROIs{i} = currECoGchannels(i).oldROI;
        end

        % Creating the ROI to average measure/freq functional connectivity map
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
            for iSeg = 1:length(currConnCell)
                currSegStruct = currConnCell{iSeg, 1};

                % Gets out the functional connectivity matrix depending on
                % the measure and frequency
                if contains(measure,'wPLI_')
                    if strcmp(freq, "delta") == 1
                        connMatrix = currSegStruct.wPLI_debias.delta;
                    elseif strcmp(freq, "theta") == 1
                        connMatrix = currSegStruct.wPLI_debias.theta;
                    elseif strcmp(freq, "alpha") == 1
                        connMatrix = currSegStruct.wPLI_debias.alpha;
                    elseif strcmp(freq, "beta") == 1
                        connMatrix = currSegStruct.wPLI_debias.beta;
                    elseif strcmp(freq, "gamma") == 1
                        connMatrix = currSegStruct.wPLI_debias.gamma;
                    end
                elseif contains(measure,'envCorr')
                    if strcmp(freq, "delta") == 1
                        connMatrix = currSegStruct.envCorr.delta;
                    elseif strcmp(freq, "theta") == 1
                        connMatrix = currSegStruct.envCorr.theta;
                    elseif strcmp(freq, "alpha") == 1
                        connMatrix = currSegStruct.envCorr.alpha;
                    elseif strcmp(freq, "beta") == 1
                        connMatrix = currSegStruct.envCorr.beta;
                    elseif strcmp(freq, "gamma") == 1
                        connMatrix = currSegStruct.envCorr.gamma;
                    elseif strcmp(freq, "highGamma") == 1
                        connMatrix = currSegStruct.envCorr.highGamma;
                    end
                else
                    error("Measure not contianed in code --> update code")
                end


                connValues = connMatrix(rowIndices, :);
                %Turns it to be 1D
                connValues = connValues(:);
                allMeanValues = [allMeanValues mean(connValues, 'omitnan')];
            end

            allMeanValues = allMeanValues(:);
            ROI2Avg(currROI) = mean(allMeanValues, 'omitnan');
            
        %end of Contacts loop
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