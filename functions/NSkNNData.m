function [dataImputed dataImputedWeighted] = NSkNNData(dataMV,K)
% Function to impute missing values in a dataset using NSkNN. dataMV is 
% the dataset with the missing values (can either be with or without 
% autoscale) and K is the # of nearest neighbors to use to impute the data.
% NSkNNData does not skip neighbors that have NaN values in the same 
% location as the metabolite being imputed. Instead, it replaces these NaN
% values with the minimum value of that metabolite.

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
        %distance = pdist2(MVRowsRemoved(:,targetCol)',MVRowsRemoved(:,neighborCol)');
        distIdx(counter,:) = [targetCol distance neighborCol];
        counter = counter+1;
    end
end

% Remove rows that calculated the Euclidean distance between a sample and
% itself.
sameSample = find(distIdx(:,1)==distIdx(:,3));
distIdx(sameSample,:) = [];
distIdxSorted = sortrows(distIdx);

minValperRow = min(dataMV,[],2); % Finds minimum of each metabolite

% Implement NSkNN with minimum value replacement
dataImputed = dataMV;
dataImputedWeighted = dataMV;
for targetCol = 1:numCol
    numMV = size(rowMV{targetCol},1); % # of missing values in a column
    firstNNIdx = (targetCol-1)*(numCol-1)+1; % Column index of the first nearest neighbor of target sample (i) in distIdxSorted
    for MVidx = 1:numMV % For each missing value in the target sample...
        tempDataMV = dataMV;
        NN = distIdxSorted(firstNNIdx:firstNNIdx+K-1,3); % Column #s of the k nearest neighbors
        DistanceNN = distIdxSorted(firstNNIdx:firstNNIdx+K-1,2); % Distances of k nearest neighbors
        idxNaNinCol = find(isnan(tempDataMV(rowMV{targetCol}(MVidx),NN))); % Finds missing values that are the same metabolite as the target metabolite to be imputed
        if isempty(idxNaNinCol)~=1 % If NaN values found...
            % If there are NaN values in the nearest neighbor
            % metabolite that is the same as the target metabolite to
            % be imputed, replace with min value of the target
            % metabolite

            tempDataMV(rowMV{targetCol}(MVidx),NN(idxNaNinCol)) = minValperRow(rowMV{targetCol}(MVidx));           
        end
        % Imputed data is weighted by the inverse of the distance
        WeightMultiplier = (1./DistanceNN')/sum(1./DistanceNN);
        dataImputedWeighted(rowMV{targetCol}(MVidx),targetCol) = sum(tempDataMV(rowMV{targetCol}(MVidx),NN).*WeightMultiplier);
        % Not weighted
        dataImputed(rowMV{targetCol}(MVidx),targetCol) = mean(tempDataMV(rowMV{targetCol}(MVidx),NN));
    end
end
