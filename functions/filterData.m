function [filteredMV filteredNoMV rowToDelete colToDelete] = filterData(dataMV,rawData)
% This function removes any samples and metabolites that have >70% (by
% default) of missing values.

numRow = size(dataMV,1);
numCol = size(dataMV,2);
percentThreshold = 70; % Percent of missing values allowed to be in row/column

% Find rows that have > percentThreshold missing
percentMVinRow = sum(isnan(dataMV),2)/numCol;
rowToDelete = find(percentMVinRow>percentThreshold/100);

% Find columns that have > percentThreshold missing
percentMVinCol = sum(isnan(dataMV),1)/numRow;
colToDelete = find(percentMVinCol>percentThreshold/100);

filteredMV = dataMV;
filteredNoMV = rawData;

% Delete identified rows
filteredMV(rowToDelete,:) = [];
filteredNoMV(rowToDelete,:) = [];

% Delete identified columns
filteredMV(:,colToDelete) = [];
filteredNoMV(:,colToDelete) = [];
