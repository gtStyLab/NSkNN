function dataMV = removeDataMCAR(rawData,percentMV)
% Function to create missing values in the original dataset (rawData)
% based on the percent missing values (percentMV) given. Values are removed
% randomly to simulate MCAR values. Output is a amatrix with missing values
% (dataMV).

numVal = numel(rawData); % # of total values in raw data
numValRemove = round(percentMV/100*numVal); % # of values to remove

% Set 0 values to small number to avoid dividing by 0 in future steps
rawData(rawData==0) = 0.00001;

% Randomly select values to remove and set them as NaN
removeIdx = randperm(numVal,numValRemove);
dataMV = rawData;
dataMV(removeIdx) = NaN;