function dataMV = removeDataMNART(rawData,percentMV,percentMNAR)
% Function to create missing values in the original data set (rawData)
% using a MNAR-T model. percentMV is the total % missingness and
% percentMNAR is the percentage of MNAR values out of all values. Output is 
% a matrix with missing values (dataMV).

% Adjust percentMNAR to account for half of values randomly removed in upper
% 20% of abundThresh.
abundThresh = percentMNAR/100/0.9; % Threshold to create MNAR values

% Set 0 values to small number to avoid dividing by 0 in future steps
rawData(rawData==0) = 0.00001;

percentRemoveAboveThresh = percentMV/100 - abundThresh*0.9;

[Abund sortAbundidx] = sort(rawData,2);
numSampBelowThresh = round(abundThresh*size(Abund,2)); % # of samples that can have MNAR values
lowAbundmet = Abund(:,1:numSampBelowThresh); % Values that can be MNAR
highAbundmet = Abund(:,numSampBelowThresh+1:end); % Values that can be MCAR

alwaysMV = round(0.8*size(lowAbundmet,2)); % # of samples that will always be MNAR
lowAbundmet(:,1:alwaysMV) = NaN; % Values that are always MNAR

% Randomly selects 50% of remaining possible LOD samples
if alwaysMV < size(lowAbundmet,2)
    lowAbundmetRandomMV = lowAbundmet(:,alwaysMV+1:end);
    lowAbundmet(:,alwaysMV+1:end) = [];

    numVal = numel(lowAbundmetRandomMV);
    numValRemove = round(0.5*numVal);
    removeIdx = randperm(numVal,numValRemove);
    lowAbundmetRandomMV(removeIdx) = NaN;

    lowAbundmet = [lowAbundmet lowAbundmetRandomMV];
end

% Start to remove MCAR values
numVal = numel(highAbundmet);
numValRemove = round(percentRemoveAboveThresh*numVal);

removeIdx = randperm(numVal,numValRemove);
highAbundmet(removeIdx) = NaN;

% Re-orders everything into original form
dataMV = [lowAbundmet highAbundmet];
for i = 1:size(dataMV,1)
    temp = sortrows([sortAbundidx(i,:)' dataMV(i,:)'])';
    dataMV(i,:) = temp(2,:);
end