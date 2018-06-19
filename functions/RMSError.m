function [RMSE,NormalizedRMS] = RMSError(imputedData,filteredMV,filteredNoMV)
% This function calculates the root mean square error (RMSE) and normalized root
% mean square error (NRMSE)

totNumMV = sum(sum(isnan(filteredMV)));
MeanOfMet = mean(filteredNoMV,2)*ones(1,size(filteredNoMV,2));
RMSE = sqrt(sum(sum((imputedData-filteredNoMV).^2))/totNumMV);
NormalizedRMS = sqrt(sum(sum(((imputedData-filteredNoMV)./MeanOfMet).^2))/totNumMV);