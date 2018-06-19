function [scaledMV avgMV stddevMV] = preprocessData(dataMV)
% preprocessData is a function that takes the missing value dataset and 
% autoscales the data by subtracting by the mean and dividing by the 
% standard deviation of each metabolite.

avgMV = nanmean(dataMV,2);
stddevMV = nanstd(dataMV,0,2);

scaledMV = (dataMV - avgMV*ones(1,size(dataMV,2)))./(stddevMV*ones(1,size(dataMV,2)));


