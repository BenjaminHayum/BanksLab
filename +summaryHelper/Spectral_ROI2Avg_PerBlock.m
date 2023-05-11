function [roi2avg] = Spectral_ROI2Avg_PerBlock(dataTable, band)
%% Documentation
%
%   Create a ROI to average specified frequency band spectral power map
%       For each ROI, get all spectral power values in the specified 
%       frequency band in a segment. Then, average across that segment. 
%       Then, average across all segments in a block. 
%       Finally, store this value in the map 
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

%% Converting it to a table to output
roi2avg = cell2table(blockCells, 'VariableNames', {'patientID', 'mapCondition', 'ROI2Avg'});

end