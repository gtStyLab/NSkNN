%% NS-kNN vs KNN-TN (MATLAB portion) (Figure 6, S7, S8)

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

K = 10;

for dataSet = 1:3
    for percentMV = [9 15 30]
        for MNAR = [1/3 2/3]
            mnarPer = MNAR*percentMV;
            
            fprintf('Working on %s, %d%% missing, %d/3 MNAR...\n',dataFile{dataSet},percentMV,MNAR*3);
            
            for iter = 1:100
                tic
                fileName = sprintf('%s_PercentMV-%02d_PercentMNAR-%02d_rep-%03d.csv',dataFile{dataSet},percentMV,mnarPer,iter);
                load(dataFile{dataSet});
                dataMV = load(fileName);

                % filterData function removes any samples that have greater than 70% of
                % missing values
                [filteredMV filteredNoMV] = filterData(dataMV,rawData);

                % preprocessData log transforms the data and normalizes the data by
                % subtracting the mean of the metabolites for each metabolite value and
                % divides by the standard deviation of the metabolite.
                [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

                % kNN on scaled data
                [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
                % NSkNN on scaled data
                [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

                imputedDataWeighted_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
                imputedDataWeighted_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

                % Calculate error
                [RMS_kNN NormalizedRMS_kNN] = RMSError(imputedDataWeighted_kNN_final,filteredMV,filteredNoMV);
                [RMS_NSkNN NormalizedRMS_NSkNN] = RMSError(imputedDataWeighted_NSkNN_final,filteredMV,filteredNoMV);

                RMS_kNN_compiled(iter,1) = mean(RMS_kNN);
                NormalizedRMS_kNN_compiled(iter,1) = mean(NormalizedRMS_kNN);
                RMS_NSkNN_compiled(iter,1) = mean(RMS_NSkNN);
                NormalizedRMS_NSkNN_compiled(iter,1) = mean(NormalizedRMS_NSkNN);
                fprintf('Iteration %d: ',iter);
                toc
            end
            saveName = sprintf('%sresults',fileName(1:end-11));
            save(saveName);
            clearvars -except dataSet dataFile percentMV MNAR K
        end
    end
end