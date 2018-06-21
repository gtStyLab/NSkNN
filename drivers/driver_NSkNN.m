% General driver for using NSkNN

clear;clc;

fileName = 'DataFile.csv';
csvData = importdata(fileName);

rawData = csvData.data;
dataMV = rawData;

K = 6;

% filterData function removes any samples or metabolites that have 
% greater than 70% of the values missing. Inputs are the missing value
% dataset (dataMV) and the original raw dataset (rawData)
[filteredMV filteredNoMV rowsRemoved colsRemoved] = filterData(dataMV,rawData);

% preprocessData normalizes the data by subtracting the mean of the 
% metabolites for each metabolite value and dividing by the standard 
% deviation of the metabolite. Input is the missing value data aftering 
% filteirng (filteredMV)
[scaledMV avgMV stddevMV] = preprocessData(filteredMV);

% NS-kNN performed on the preprocessed data (scaledMV)
[imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

% Return imputed data to original scale by post processing the data
imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);