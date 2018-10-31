clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData','UrinaryHumanData','NeuralHumanData','MicrobeData','AntibioticMouseData','TobaccoData','SMData'};

K = 6;

for dataSet = 1:9
    for percentMV = [10 30]
        for percMNAR = [25 50 75]
            for percentMVlowAbund_III = [30 40]
                fprintf('Working on %s, %d%% missing, %d%% MNAR, Threshold III = %d...\n',dataFile{dataSet},percentMV,percMNAR,percentMVlowAbund_III);

                dataName = dataFile{dataSet};
                load(dataName);
            
                for i = 1:100

                    fileName = sprintf('%s_MM_PercMV-%02d_ThreshIII-%02d_PercMNAR-%02d_rep-%03d.csv',dataFile{dataSet},percentMV,percentMVlowAbund_III,percMNAR,i);
                    dataMV = load(fileName);

                    % filterData function removes any samples that have greater than 70% of
                    % missing values
                    [filteredMV filteredNoMV] = filterData(dataMV,rawData);

                    % preprocessData log transforms the data and normalizes the data by
                    % subtracting the mean of the metabolites for each metabolite value and
                    % divides by the standard deviation of the metabolite.
                    [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

                    % kNN on scaled data (Euclidean)
                    [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
                    % kNN on scaled data (Pearson)
                    [imputedData_kNNPearson imputedDataWeighted_kNNPearson] = kNNData_PearsonCorr(scaledMV,K);
                    % NSkNN on scaled data
                    [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

                    imputedDataWeighted_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
                    imputedDataWeighted_kNNPearson_final = postprocessData(imputedDataWeighted_kNNPearson, avgMV, stddevMV);
                    imputedDataWeighted_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

                    % Calculate error
                    [RMS_kNN NormalizedRMS_kNN] = RMSError(imputedDataWeighted_kNN_final,filteredMV,filteredNoMV);
                    [RMS_kNNPearson NormalizedRMS_kNNPearson] = RMSError(imputedDataWeighted_kNNPearson_final,filteredMV,filteredNoMV);
                    [RMS_NSkNN NormalizedRMS_NSkNN] = RMSError(imputedDataWeighted_NSkNN_final,filteredMV,filteredNoMV);

                    RMS_kNN_compiled(i,1) = mean(RMS_kNN);
                    NormalizedRMS_kNN_compiled(i,1) = mean(NormalizedRMS_kNN);
                    RMS_kNNPearson_compiled(i,1) = mean(RMS_kNNPearson);
                    NormalizedRMS_kNNPearson_compiled(i,1) = mean(NormalizedRMS_kNNPearson);
                    RMS_NSkNN_compiled(i,1) = mean(RMS_NSkNN);
                    NormalizedRMS_NSkNN_compiled(i,1) = mean(NormalizedRMS_NSkNN);

                end
                compiledError = [NormalizedRMS_NSkNN_compiled NormalizedRMS_kNN_compiled NormalizedRMS_kNNPearson_compiled];
                saveName = sprintf('results\\%sMATLABresults.csv',fileName(1:end-11));
                csvwrite(saveName,compiledError);
            end
        end
    end
end



