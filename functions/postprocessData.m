function postprocessImputed = postprocessData(imputedData, avgMV, stddevMV)
% The postprocessData function returns the imputed data into the
% original scale of the raw data.


postprocessImputed = imputedData.*(stddevMV*ones(1,size(imputedData,2)))+avgMV*ones(1,size(imputedData,2));