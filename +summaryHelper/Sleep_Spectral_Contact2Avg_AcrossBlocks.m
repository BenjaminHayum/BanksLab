function [contact2avg] = Sleep_Spectral_Contact2Avg_AcrossBlocks(dataTable, band)
%% Documentation
%
%   Create a contact to spectral band to average spectral power map
%       Get all spectral values for the band in that segment.
%       Then, average across that segment. 
%       Then, average across all segments in a block. 
%       Finally, store this value in the map 
%
%   ***ONLY ONE BLOCK IN SLEEP DATASET SO THERE'S NO NEED TO AVERAGE
%   ACROSS MULTIPLE BLOCKS***
%
%   Make sure the data table includes the mapped condition names already 
%       --> make a table column called "mapCondition"
%
%   ***ONLY WORKS FOR N2/N3 and WS. DOES NOT GO OVER EVERY SLEEP STATE.
%   SPECIFICALLY DESIGN TO BE USED IN THE SPECTRAL RATIO PLOT V2's***
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

%% Gathering the Sleep Spectral Data by ROI
% Column 1 is Subject
% Column 2 is whether subject is N2 or N3
% Column 3 is N3/N2 band average spectral power Channel dictionary
% Column 4 is WS band average spectral power Channel dictionary
sleepCells = cell(length(dataTable.patientID), 4);
[patientList, ia, ib] = unique(dataTable.patientID,'stable');

dictionaryCounter = 0;
for iPatient = 1:length(patientList)
    currPatient = patientList{iPatient};
    dictionaryCounter = dictionaryCounter + 1;

    patientSpecCells = dataTable.specAnalysis{iPatient};

    % Listing the states for this patient
    patientStates = dataTable.states(iPatient);
    patientStates = patientStates{:};
    allStates = unique(patientStates);

    N2Boolean = 0;
    if sum(strcmp(allStates, "N3")) == 0
        N2Boolean = 1;
    end
    for iCondition = 1:length(allStates)
        currState = allStates{iCondition};
        if strcmp(currState, "WS") == 1
            currState = "W";
        end
        
        % Initializing the Maps (to be used in DictionaryCells)
        conditionMap = containers.Map();

        % Getting the Channels pre-organized
        currECoGchannels = dataTable.ECoGchannels{iPatient};
        Channels = cell(1, length(currECoGchannels));
        for i = 1:length(currECoGchannels)
            Channels{i} = num2str(currECoGchannels(i).contNum);
        end

        % Loop through all unique Chanels and add to dictionary
        for iChannel = 1:length(Channels)
            % Get row indices
            currChannel = Channels{iChannel};
            rowIndex = iChannel;
            
            %Now, index into spectralCells and append to allValues
            allMeanValues = [];
            currSpecCell = patientSpecCells;
            specNames = fieldnames(currSpecCell);
            for iSegStruct = 1:size(specNames)
                if contains(specNames{iSegStruct}, currState) == 1
                    currSegStruct = currSpecCell.(specNames{iSegStruct});
                
                    % Based off of column indices in "freq" that correspond to the
                    % frequencies recorded in "powspctrm"
                    columnIndices = [];
                    for i = 1:length(currSegStruct.freq)
                        if (currSegStruct.freq(i) >= bandFreq(1)) && (currSegStruct.freq(i) <= bandFreq(2))
                            columnIndices = [columnIndices i];
                        end
                    end
    
                    specValues = currSegStruct.powspctrm(rowIndex, columnIndices);
                    allMeanValues = [allMeanValues mean(specValues(:), 'omitnan')];
                end
            end
            %making a map of each ROI and its corresponding mean spectral
            %power
            conditionMap(currChannel) = mean(allMeanValues(:), 'omitnan');
            
        %end of allROIs loop
        end
       
        % Putting it all into the cells

        % Current patient
        sleepCells{dictionaryCounter, 1} = currPatient;
        % N2 vs N3 status
        if N2Boolean == 1
            sleepCells{dictionaryCounter, 2} = "N2";
        else
            sleepCells{dictionaryCounter, 2} = "N3";
        end
        % Entering the dictionaries
        if strcmp(currState, "W") == 1
            sleepCells{dictionaryCounter, 4} = conditionMap;
        elseif strcmp(currState, "N3") == 1
            sleepCells{dictionaryCounter, 3} = conditionMap;
        elseif strcmp(currState, "N2") == 1 && N2Boolean == 1
            sleepCells{dictionaryCounter, 3} = conditionMap;
        end

    %End of states loop
    end

    %End of patients loop
end
%% Output

contact2avg = cell2table(sleepCells, 'VariableNames', {'patientID', 'state', 'N3N2contact2avg', 'WScontact2avg'});

end