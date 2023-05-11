function [roi2avg] = FuncConn_ROI2Avg_PerBlock(dataTable, measure, freq)
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

%% Converting it to a table to output
roi2avg = cell2table(blockCells, 'VariableNames', {'patientID', 'mapCondition', 'ROI2Avg'});

end