function [contact2avg] = FuncConn_Contact2Avg_AcrossBlocks(dataTable, measure, freq)
%% Documentation
%
%   Create a contact to average specified measure and frequency functional 
%       connectivity map
%
%   For each contact, get all functional connectivity values in a segment. 
%       Then, average across that segment. 
%       Then, average across all segments in a block. 
%       Finally, store this value in the map 
%
%   ***THEN, average across all subjects with multiple blocks of the same
%       condition***
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
% 3rd column is Contact/Channel to Average Measure/Freq Functional Connectivity Map
% 4th column is Contact/Channel to ROI map

[patientList, ia, ib] = unique(dataTable.patientID,'stable');

blockCells = cell(length(dataTable.patientID), 4);
blockCounter = 0;
for iPatient = 1:length(patientList)
    patientConditions = dataTable.mapCondition(ib==iPatient);
    for iCondition = 1:length(patientConditions)
        blockCounter = blockCounter + 1;
        currPatient = patientList{iPatient};
        currCondition = dataTable.mapCondition{blockCounter};
        currECoGchannels = dataTable.ECoGchannels{blockCounter};
        currConnCell = dataTable.(measure){blockCounter};

        % Getting list of all Contacts and mapping their corresponding ROI
        Contact2ROI = containers.Map();
        Contacts = cell(1, length(currECoGchannels));
        for i = 1:length(currECoGchannels)
            currContact = num2str(currECoGchannels(i).contNum);
            currROI = currECoGchannels(i).oldROI;
            Contacts{i} = currContact;
            Contact2ROI(currContact) = currROI;
        end

        % Creating the contact to average measure/freq functional connectivity map
        %   Extract all values out of a segment, average across that
        %   segment, then average across all segments
        Contact2Avg = containers.Map();
        for iContact = 1:length(Contacts)
            currContact = Contacts{iContact};
            rowIndex = iContact;
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


                connValues = connMatrix(rowIndex, :);
                %Turns it to be 1D
                connValues = connValues(:);
                allMeanValues = [allMeanValues mean(connValues, 'omitnan')];
            end

            allMeanValues = allMeanValues(:);
            Contact2Avg(currContact) = mean(allMeanValues, 'omitnan');
            
        %end of Contacts loop
        end

        blockCells{blockCounter, 1} = currPatient;
        blockCells{blockCounter, 2} = currCondition;
        blockCells{blockCounter, 3} = Contact2Avg;
        blockCells{blockCounter, 4} = Contact2ROI;
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
% 3rd column is contact to avg across blocks map
% 4th column is contact to ROI map

comboCells = cell(uniqueCounter, 4);
comboCounter = 0;
patientList = unique(dataTable.patientID,'stable');
for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    currConditions = patient2conditions(currPatient);
    for iCondition = 1:length(currConditions)
        currCondition = currConditions{iCondition};

        % To store the array of all averages
        arrayChannelMap = containers.Map();
        % To store the averages of the averages across subjects
        outputChannelMap = containers.Map();

        % Find all cells with the current patient & condition combination
        for iCell = 1:length(blockCells)
            if strcmp(currPatient, blockCells{iCell, 1}) == 1 && strcmp(currCondition, blockCells{iCell, 2}) == 1
                currMap = blockCells{iCell, 3};
                allChannels = keys(currMap);
                for iChannel = 1:length(allChannels)
                    currChannel = allChannels{iChannel};
                    currAvg = currMap(currChannel);

                    % If the contact is already in arrayContactMap, add
                    % this average value to the array stored there
                    %
                    % Else, create the start of a new average array for 
                    % this contact
                    if any(strcmp(currChannel, keys(arrayChannelMap))) == 1
                        currAvgArray = arrayChannelMap(currChannel);
                        currAvgArray = [currAvgArray currAvg];
                        arrayChannelMap(currChannel) = currAvgArray;
                    else
                        arrayChannelMap(currChannel) = currAvg;
                    end
                end

                Channel2ROI = blockCells{iCell, 4};
            end
        end
        %End of blockCells loop
        
        % Calculating averages and assigning to plotCells
        for iChannel = 1:length(allChannels)
            currChannel = allChannels{iChannel};
            currAvgArray = arrayChannelMap(currChannel);
            outputChannelMap(currChannel) = mean(currAvgArray, 'omitnan');
        end
        
        comboCounter = comboCounter + 1;
        comboCells{comboCounter, 1} = currPatient;
        comboCells{comboCounter, 2} = currCondition;
        comboCells{comboCounter, 3} = outputChannelMap;
        comboCells{comboCounter, 4} = Channel2ROI;

    end
    %End of conditions loop
end
% End of patients loop


%% Converting it to a table to output
contact2avg = cell2table(comboCells, 'VariableNames', {'patientID', 'mapCondition', 'Contact2Avg', 'Contact2ROI'});

end