function [imputedData imputedDataWeighted] = kNNData(dataMV,K)
% Function to impute missing values in a dataset using NSkNN. dataMV is 
% the dataset with the missing values (can either be with or without 
% autoscale) and K is the # of nearest neighbors to use to impute the data.
% kNNData skips neighbors that have NaN values in the same location as the 
% metabolite being imputed and move on to the next nearest neighbor.

numCol = size(dataMV,2);
for col = 1:numCol
    rowMV{col} = find(isnan(dataMV(:,col))); % Finds the row # of every missing value in each column
end

counter = 1;
for targetCol = 1:numCol % i is the target sample
    for neighborCol = 1:numCol % Calculate the Euclidean distances between the target sample (i) and the other samples (j)
        MVRowsRemoved = dataMV;
        rowsToRemove = union(rowMV{targetCol},rowMV{neighborCol}); % Ignore NaNs when calculating distances
        MVRowsRemoved(rowsToRemove,:) = []; % Remove rows in target sample that have missing values
        numMetInCalc = size(MVRowsRemoved,1); % # of metabolites used in calculation of metabolites in order to weight distance
        % Divide by numMetInCalc to avoid scenarios where distances 
        % calculated with only a few metabolites are weighted more heavily 
        % over distances that are close to the same distance, but used more 
        % metabolites in the calculation.
        distance = pdist2(MVRowsRemoved(:,targetCol)',MVRowsRemoved(:,neighborCol)')/sqrt(numMetInCalc);         
        distIdx(counter,:) = [targetCol distance neighborCol];
        counter = counter+1;
    end
end

% Remove rows that calculated the Euclidean distance between a sample and
% itself.
sameSample = find(distIdx(:,1)==distIdx(:,3));
distIdx(sameSample,:) = [];
distIdxSorted = sortrows(distIdx);


imputedData = dataMV;
imputedDataWeighted = dataMV;

for targetCol = 1:numCol
    numMV = size(rowMV{targetCol},1); % # of missing values in the target sample
    firstNNIdx = (targetCol-1)*(numCol-1)+1; % Column index of the first nearest neighbor of target sample (targetCol) in distIdxSorted
    NNCleared = zeros(numMV,1); % # of NN that are cleared of matching missing values with the target metabolite
    while sum(NNCleared) < numMV
        % Keep looping if there are missing values in the same location for
        % the value being imputed and the nearest neighbors
        for MVidx = 1:numMV % For each missing value in the target sample...
            % Find NN that have a missing value for the same metabolite
            % that is being imputed.
            Ktemp = K;
            NN = distIdxSorted(firstNNIdx:firstNNIdx+Ktemp-1,3); % Column # of the k nearest neighbors
            DistanceNN = distIdxSorted(firstNNIdx:firstNNIdx+Ktemp-1,2); % Distance of nearest neighbors
            while NNCleared(MVidx)==0
                idxNaNinCol = find(isnan(dataMV(rowMV{targetCol}(MVidx),NN)));
                NumNNtoReplace = size(idxNaNinCol,2); % # of NN to replace
                if isempty(idxNaNinCol)~=1 % If NaN values found...
                    NN(idxNaNinCol) = []; % Remove NN that have NaN at target metabolite            
                    DistanceNN(idxNaNinCol) = [];
                    NN = [NN; distIdxSorted(firstNNIdx+Ktemp:firstNNIdx+Ktemp+NumNNtoReplace-1,3)]; % Add new NNs
                    DistanceNN = [DistanceNN; distIdxSorted(firstNNIdx+Ktemp:firstNNIdx+Ktemp+NumNNtoReplace-1,2)];
                    Ktemp = Ktemp + NumNNtoReplace;
                else
                    % If no missing values found in columns that match the
                    % missing metabolite value being imputed, the nearest
                    % neighbor is cleared.
                    NNCleared(MVidx)=1;
                end
            end
            % Imputed data is weighted by the inverse of the distance
            WeightMultiplier = (1./DistanceNN')/sum(1./DistanceNN);
            imputedDataWeighted(rowMV{targetCol}(MVidx),targetCol) = sum(dataMV(rowMV{targetCol}(MVidx),NN).*WeightMultiplier);
            % Imputed data is not weighted
            imputedData(rowMV{targetCol}(MVidx),targetCol) = mean(dataMV(rowMV{targetCol}(MVidx),NN));
        end
    end
end
