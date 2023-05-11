function [currMin,currMax] = squareLogAxisLimits(allXPlotted,allYPlotted)
%% Documentation
%
% Getting the axes limits for a square log plot
%
% For square log axes limits, you need a minimum > 0 
% Maximum should naturally be > 0
%

%% Getting the axes limits

firstMin = min(allXPlotted);
firstMax = max(allXPlotted);
secondMin = min(allYPlotted);
secondMax = max(allYPlotted);
if secondMin < firstMin
    currMin = secondMin;
else
    currMin = firstMin;
end

% To make the plots work with the log scale, we can't set min to <= 0
% So, get the smallest > 0 value:
if currMin < 0
    sortedFirst = sort(allFirstPlotted, 'ascend');
    for iVal = 1:length(sortedFirst)
        if sortedFirst(iVal) > 0
            tempFirstMin = sortedFirst(iVal);
            break;
        end
    end
    sortedSecond = sort(allSecondPlotted, 'ascend');
    for iVal = 1:length(allSecondPlotted)
        if sortedSecond(iVal) > 0
            tempSecondMin = sortedSecond(iVal);
            break;
        end
    end
    if tempFirstMin < tempSecondMin
        currMin = tempFirstMin;
    elseif tempSecondMin <= tempFirstMin
        currMin = tempSecondMin;
    end
end


if secondMax > firstMax
    currMax = secondMax;
else
    currMax = firstMax;
end

end