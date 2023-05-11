function [contact2avg] = FuncConn_Contact2Avg_PerBlock(dataTable, measure, freq)
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
%   Have accomodating contact to ROI map for later color coded plots   
%
%   Make sure the data table includes the mapped condition names already 
%       --> make a table column called "mapCondition"
%


%% Creating the Contact2Average
% FOR EVERY SINGLE BLOCK IN THE DATA TABLE:
% 1st column is patient
% 2nd column is condition
% 3rd column is Contact/Channel to Average Band Spectral Power Map
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

%% Converting it to a table to output
contact2avg = cell2table(blockCells, 'VariableNames', {'patientID', 'mapCondition', 'Contact2Avg', 'Contact2ROI'});

end