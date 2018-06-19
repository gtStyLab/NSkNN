function dataMV = removeDataMM(rawData,percentMV,percentBelowThresh_I,LowAbundThresh_II,percentMVlowAbund_III)
% Function to create missing values in the original data set (rawData)
% using a Mixed Missingness (MM) model. percentMV is the total %
% missingness. Descriptions of percentBelowThresh_I, LowAbundThresh_II, and
% percentMVlowAbund_III can be found in the Methods section. Output is a 
% matrix with missing values (dataMV).

%rng('shuffle');
percentMVlowAbund_III = percentMVlowAbund_III/0.9;
percentMVmedAbund = 0.5*percentMVlowAbund_III;

TotalNumValRemove = round(percentMV/100*numel(rawData)); % # of values to remove

% Set 0 values to small number to avoid dividing by 0 in future steps
rawData(rawData==0) = 0.00001;

% Finds average abundance across all metabolites and sorts from
% least amount to greatest amount
[avgMetabolite sortAvgIdx] = sort(nanmean(rawData,2),'ascend');

numMetBelowThresh = round(percentBelowThresh_I/100*size(avgMetabolite,1)); % # of metabolites that have an avg abundance below the threshold

belowThreshMet = rawData(sortAvgIdx(1:numMetBelowThresh),:); % Create matrix of values that are below threshold
belowThreshIdx = sortAvgIdx(1:numMetBelowThresh); % Row indices from the raw data for the values in belowThreshMet
aboveThreshMet = rawData(sortAvgIdx(numMetBelowThresh+1:end),:); % Create matrix of values that are above threshold
aboveThreshIdx = sortAvgIdx(numMetBelowThresh+1:end); % Row indices from the raw data for the values in aboveThreshMet

% Split metabolites below initial threshold further into low abundance and medium abundance
numMetLowAbund = round(LowAbundThresh_II/100*size(belowThreshMet,1));
lowAbundMet = rawData(sortAvgIdx(1:numMetLowAbund),:);
medAbundmet = rawData(sortAvgIdx(numMetLowAbund+1:size(belowThreshMet,1)),:);


% # of samples for each metabolite that belong to the low abundance group of metabolites that
% will be missing.
numSampBelowAbundThreshLowAbund = round(percentMVlowAbund_III/100*size(lowAbundMet,2));
[sortLowAbundmet sortAbundLowAbundidx] = sort(lowAbundMet,2); % Sort matrix with low mean abundance metabolites from low to high abundance values for each metabolite.
LowLowAbund = sortLowAbundmet(:,1:numSampBelowAbundThreshLowAbund); % Create matrix with low abundance values for each metabolite with a low mean abundance.
HighLowAbund = sortLowAbundmet(:,numSampBelowAbundThreshLowAbund+1:end); % Create a matrix with high abundance values for each metabolite with a low mean abundance.

alwaysMV = round(0.8*size(LowLowAbund,2));
LowLowAbund(:,1:alwaysMV) = NaN;
if alwaysMV < size(LowLowAbund,2)
    LowLowAbundRandomMV = LowLowAbund(:,alwaysMV+1:end);
    LowLowAbund(:,alwaysMV+1:end) = [];
    
    numVal = numel(LowLowAbundRandomMV);
    numValRemove = round(0.5*numVal);
    removeIdx = randperm(numVal,numValRemove);
    LowLowAbundRandomMV(removeIdx) = NaN;

    LowLowAbund = [LowLowAbund LowLowAbundRandomMV];
end


% # of samples for each metabolite that belong to the medium abundance group of metabolites that
% will be missing.
numSampBelowAbundThreshMedAbund = round(percentMVmedAbund/100*size(medAbundmet,2));
[sortmedAbundmet sortAbundMedAbundidx] = sort(medAbundmet,2); % Sort matrix with high mean abundance metabolites from low to high abundance values for each metabolite.
LowMedAbund = sortmedAbundmet(:,1:numSampBelowAbundThreshMedAbund); % Create matrix with low abundance values for each metabolite with a high mean abundance.
HighMedAbund = sortmedAbundmet(:,numSampBelowAbundThreshMedAbund+1:end); % Create a matrix with high abundance values for each metabolite with a high mean abundance.

alwaysMV = round(0.8*size(LowMedAbund,2));
LowMedAbund(:,1:alwaysMV) = NaN;
if alwaysMV < size(LowLowAbund,2)
    LowMedAbundRandomMV = LowMedAbund(:,alwaysMV+1:end);
    LowMedAbund(:,alwaysMV+1:end) = [];
    
    numVal = numel(LowMedAbundRandomMV);
    numValRemove = round(0.5*numVal);
    removeIdx = randperm(numVal,numValRemove);
    LowMedAbundRandomMV(removeIdx) = NaN;

    LowMedAbund = [LowMedAbund LowMedAbundRandomMV];
end

% Add MCAR values
numMVbelowThreshLOD = sum(sum(isnan(LowMedAbund))) + sum(sum(isnan(LowLowAbund)));
numMVrandom = TotalNumValRemove - numMVbelowThreshLOD;
if numMVrandom < 0
    numMVrandom = 0;
end
MVrandomAvailable = numel(rawData)-numel(LowMedAbund)-numel(LowLowAbund);

aboveThreshWithMV = aboveThreshMet;
HighLowAbundWithMV = HighLowAbund;
HighMedAbundWithMV = HighMedAbund;

size_aboveThreshWithMV = size(aboveThreshWithMV);
size_HighLowAbundWithMV = size(HighLowAbundWithMV);
size_HighMedAbundWithMV = size(HighMedAbundWithMV);

remainingValuesVector = [aboveThreshWithMV(:); HighLowAbundWithMV(:); HighMedAbundWithMV(:)];
removeIdx = randperm(length(remainingValuesVector),numMVrandom);
remainingValuesVector(removeIdx) = NaN;

aboveThreshWithMV = remainingValuesVector(1:numel(aboveThreshWithMV));
HighLowAbundWithMV = remainingValuesVector(numel(aboveThreshWithMV)+1:numel(aboveThreshWithMV)+numel(HighLowAbundWithMV));
HighMedAbundWithMV = remainingValuesVector(numel(aboveThreshWithMV)+numel(HighLowAbundWithMV)+1:end);

aboveThreshWithMV = reshape(aboveThreshWithMV,size_aboveThreshWithMV);
HighLowAbundWithMV = reshape(HighLowAbundWithMV,size_HighLowAbundWithMV);
HighMedAbundWithMV = reshape(HighMedAbundWithMV,size_HighMedAbundWithMV);


% Recreate original matrix with missing values and use sortrows to 
% reorder values to their original indices.
MedAbundWithMV = [LowMedAbund HighMedAbundWithMV];
for i = 1:size(MedAbundWithMV,1)
    temp = sortrows([sortAbundMedAbundidx(i,:)' MedAbundWithMV(i,:)'])';
    MedAbundWithMV(i,:) = temp(2,:);
end

LowAbundWithMV = [LowLowAbund HighLowAbundWithMV];
for i = 1:size(LowAbundWithMV,1)
    temp = sortrows([sortAbundLowAbundidx(i,:)' LowAbundWithMV(i,:)'])';
    LowAbundWithMV(i,:) = temp(2,:);
end

belowThreshWithMV = [LowAbundWithMV; MedAbundWithMV];

dataMVtemp = [belowThreshWithMV; aboveThreshWithMV];
dataMVtemp = sortrows([sortAvgIdx dataMVtemp]);
dataMV = dataMVtemp(:,2:end);