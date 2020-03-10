function dataMV = removeDataMNAR_LODparam(rawData,percentMV,LOD)
% Function to create missing values in the original data set (rawData)
% based on the percent missing values (percentMV) given. Missing values are
% created to simulate MNAR values by creating a threshold where a certain
% percent (set by the LOD input parameter) of lowest abundance values are 
% removed in each metabolite to simulate values below the limit of 
% detection (LOD) and half of the remaining (1-LOD)% are removed at random 
% to simulate missing values close to the LOD. Output is a matrix with 
% missing values (dataMV). This function is used for Figure S1 of the 
% Supplementary. In removeDataMNAR.m, the LOD input parameter is held at 
% 80%.

% Adjust percentMV to account for half of values randomly removed in upper
% portion (based on LOD parameter) of abundThresh.
LOD = LOD/100;
percDiv = (1-LOD)/2+LOD;
percentMV = percentMV/percDiv;

% Set 0 values to small number to avoid dividing by 0 in future steps
rawData(rawData==0) = 0.00001;

abundThresh = percentMV/100; % Threshold for MNAR values

[Abund sortAbundIdx] = sort(rawData,2);
numSampBelowThresh = round(abundThresh*size(Abund,2));
lowAbundMet = Abund(:,1:numSampBelowThresh);
highAbundMet = Abund(:,numSampBelowThresh+1:end);

% Remove lowest 80% values below threshold
alwaysMV = round(LOD*size(lowAbundMet,2));
lowAbundMet(:,1:alwaysMV) = NaN;

% Randomly remove half of upper 20% values below threshold
if alwaysMV < size(lowAbundMet,2)
    lowAbundMetRandomMV = lowAbundMet(:,alwaysMV+1:end);
    lowAbundMet(:,alwaysMV+1:end) = [];

    numVal = numel(lowAbundMetRandomMV);
    numValRemove = round(0.5*numVal);
    removeIdx = randperm(numVal,numValRemove);
    lowAbundMetRandomMV(removeIdx) = NaN;

    lowAbundMet = [lowAbundMet lowAbundMetRandomMV];
end

% Reorder values
dataMV = [lowAbundMet highAbundMet];
for i = 1:size(dataMV,1)
    temp = sortrows([sortAbundIdx(i,:)' dataMV(i,:)'])';
    dataMV(i,:) = temp(2,:);
end