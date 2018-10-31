%% MNAR testing LOD parameter (Figure S1)

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

K = 6;

for dataSet = 1:3
    for LOD = [20 40 60 80] % Test MNAR LOD thresholds of 20%, 40%, and 60% 80%
        
        % Load data
        fileName = dataFile{dataSet};
        load(fileName);
        fprintf('Working on %s, LOD %d%%...\n',fileName,LOD);
        
        for iter = 1:100
            tic
            for percentMV = 0:30
               
                % removeDataMNAR_LODparam function generates MNAR missing
                % values based on the total % missingness (percentMV) chosen
                % and the LOD value
                dataMV = removeDataMNAR_LODparam(rawData,percentMV,LOD);

                % filterData function removes any samples that have greater
                % than 70% of missing values
                [filteredMV filteredNoMV] = filterData(dataMV,rawData);

                % Preprocess the data by autoscaling
                [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

                % Impute with kNN
                [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
                % Impute with NSkNN
                [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

                % Postprocess the data to return data to original scale
                imputedData_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
                imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

                % Calculate error
                [RMSE_kNN(iter,percentMV+1) NormalizedRMS_kNN(iter,percentMV+1)] = RMSError(imputedData_kNN_final,filteredMV,filteredNoMV);
                [RMSE_NSkNN(iter,percentMV+1) NormalizedRMS_NSkNN(iter,percentMV+1)] = RMSError(imputedData_NSkNN_final,filteredMV,filteredNoMV);
            end
            toc
        end
        save(sprintf('%s_MNAR_%02d',fileName,LOD));
        clearvars -except dataSet dataFile K LOD
    end
end

%% MNAR-Titrate (MNAR-T) comparison of kNN, NSkNN, NSkNN_HM, NSkNN_Zero (Figure S2)

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

percentMV = 30;
K = 6;

for dataSet = 1:3

    % Load data
    fileName = dataFile{dataSet};
    load(fileName);
    fprintf('Working on %s, %d%% missing...\n',fileName,percentMV);

    for iter = 1:100
        tic
        for percentMNAR = 0:percentMV
            
            % removeDataMNART removes data using the MNAR-T modeling method
            % with inputs of the raw data (rawData), total percent
            % missingness (percentMV), and percent of all values that are
            % MNAR (percentMNAR)
            dataMV = removeDataMNART(rawData,percentMV,percentMNAR);

            % filterData function removes any samples or metabolites that have 
            % greater than 70% of the values missing. Inputs are the missing value
            % dataset (dataMV) and the original raw dataset (rawData)
            [filteredMV filteredNoMV] = filterData(dataMV,rawData);

            % preprocessData normalizes the data by subtracting the mean of the 
            % metabolites for each metabolite value and dividing by the standard 
            % deviation of the metabolite. Input is the missing value data aftering 
            % filteirng (filteredMV)
            [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

            % kNN performed on the preprocessed data (scaledMV)
            [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
            % NSkNN with minimum value replacement performed on the preprocessed data (scaledMV)
            [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);
            % NSkNN with half minimum value replacement performed on the preprocessed data (scaledMV)
            [imputedData_NSkNN_HM imputedDataWeighted_NSkNN_HM] = NSkNNData_HM(scaledMV,K,avgMV,stddevMV);
            % NSkNN with zero replacmeent performed on the preprocessed data (scaledMV)
            [imputedData_NSkNN_zero imputedDataWeighted_NSkNN_zero] = NSkNNData_Zero(scaledMV,K,avgMV,stddevMV);

            % Return imputed data to original scale by post processing the data
            imputedData_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
            imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);
            imputedData_NSkNN_HM_final = postprocessData(imputedDataWeighted_NSkNN_HM, avgMV, stddevMV);
            imputedData_NSkNN_zero_final = postprocessData(imputedDataWeighted_NSkNN_zero, avgMV, stddevMV);

            % Calculate the normalized root mean square error for each
            % imputation method
            [RMSE_kNN(iter,percentMNAR+1) NormalizedRMS_kNN(iter,percentMNAR+1)] = RMSError(imputedData_kNN_final,filteredMV,filteredNoMV);
            [RMSE_NSkNN(iter,percentMNAR+1) NormalizedRMS_NSkNN(iter,percentMNAR+1)] = RMSError(imputedData_NSkNN_final,filteredMV,filteredNoMV);      
            [RMSE_NSkNN_HM(iter,percentMNAR+1) NormalizedRMS_NSkNN_HM(iter,percentMNAR+1)] = RMSError(imputedData_NSkNN_HM_final,filteredMV,filteredNoMV);      
            [RMSE_NSkNN_zero(iter,percentMNAR+1) NormalizedRMS_NSkNN_zero(iter,percentMNAR+1)] = RMSError(imputedData_NSkNN_zero_final,filteredMV,filteredNoMV);      
        end
        fprintf('Iteration %d: ',iter);
        toc
    end
    save(sprintf('%s_MNART_MV%02d_compareNSkNNReplacement',fileName,percentMV));
    clearvars -except dataSet dataFile percentMV K
end

%% MCAR Best K value (Figure S3)

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

MaxK = 10; % Maximum K value to test
percentMV = 30; % Set percent of missing values to be taken out of raw data. Code breaks if >~40%

for dataSet = 1:3
    
    % Load data
    fileName = dataFile{dataSet};
    load(fileName);
    fprintf('Working on %s, %d%% missing...\n',fileName,percentMV);
    
    for iter = 1:100
        tic
        
        % removeDataMCAR function generates a missing value dataset with only 
        % MCAR values. Inputs are the raw dataset (rawData) and total 
        % missingness percentage (percentMV)
        dataMV = removeDataMCAR(rawData,percentMV);

        % filterData function removes any samples or metabolites that have 
        % greater than 70% of the values missing. Inputs are the missing value
        % dataset (dataMV) and the original raw dataset (rawData)
        [filteredMV filteredNoMV] = filterData(dataMV,rawData);

        % preprocessData normalizes the data by subtracting the mean of the 
        % metabolites for each metabolite value and dividing by the standard 
        % deviation of the metabolite. Input is the missing value data aftering 
        % filteirng (filteredMV)
        [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

        for K=1:MaxK
            % kNN performed on the preprocessed data (scaledMV)
            [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
            % NS-kNN performed on the preprocessed data (scaledMV)
            [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

            % Return imputed data to original scale by post processing the data
            imputedData_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
            imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

            % Calculate the normalized root mean square error for each
            % imputation method
            [RMS_kNN(iter,K) NormalizedRMS_kNN(iter,K)] = RMSError(imputedData_kNN_final,filteredMV,filteredNoMV);
            [RMS_NSkNN(iter,K) NormalizedRMS_NSkNN(iter,K)] = RMSError(imputedData_NSkNN_final,filteredMV,filteredNoMV);

        end
        fprintf('Iteration %d: ',iter);
        toc
    end
    save(sprintf('%s_MCAR_K',fileName));
    clearvars -except dataSet dataFile MaxK
end

%% MCAR increasing % MV (Figure 2)

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

K = 6;

for dataSet = 1:3
    
    % Load data
    fileName = dataFile{dataSet};
    load(fileName);
    fprintf('Working on %s...\n',fileName);
    
    for iter = 1:100
        tic
        for percentMV = 0:30

        % removeDataMCAR function generates a missing value dataset with only 
        % MCAR values. Inputs are the raw dataset (rawData) and total 
        % missingness percentage (percentMV)
        dataMV = removeDataMCAR(rawData,percentMV);

        % filterData function removes any samples or metabolites that have 
        % greater than 70% of the values missing. Inputs are the missing value
        % dataset (dataMV) and the original raw dataset (rawData)
        [filteredMV filteredNoMV] = filterData(dataMV,rawData);

        % preprocessData normalizes the data by subtracting the mean of the 
        % metabolites for each metabolite value and dividing by the standard 
        % deviation of the metabolite. Input is the missing value data aftering 
        % filteirng (filteredMV)
        [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

        % kNN performed on the preprocessed data (scaledMV)
        [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
        % NS-kNN performed on the preprocessed data (scaledMV)
        [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

        % Return imputed data to original scale by post processing the data
        imputedData_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
        imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

        % Calculate the normalized root mean square error for each
        % imputation method
        [RMSE_kNN(iter,percentMV+1) NormalizedRMS_kNN(iter,percentMV+1)] = RMSError(imputedData_kNN_final,filteredMV,filteredNoMV);
        [RMSE_NSkNN(iter,percentMV+1) NormalizedRMS_NSkNN(iter,percentMV+1)] = RMSError(imputedData_NSkNN_final,filteredMV,filteredNoMV);

        end
        fprintf('Iteration %d: ',iter);
        toc
    end
    save(sprintf('%s_MCAR',fileName));
    clearvars -except dataSet dataFile K
end


%% MNAR (Figure 3)

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

K = 6;

for dataSet = 1:3

    % Load data
    fileName = dataFile{dataSet};
    load(fileName);
    fprintf('Working on %s...\n',fileName);

    for iter = 1:100
        tic
        for percentMV = 0:30
            
            % removeDataMNAR generates missing values using the MNAR
            % modeling method. Inputs are the raw data (rawData) and the
            % total percent missingness (percentMV)
            dataMV = removeDataMNAR(rawData,percentMV);

            % filterData function removes any samples or metabolites that have 
            % greater than 70% of the values missing. Inputs are the missing value
            % dataset (dataMV) and the original raw dataset (rawData)
            [filteredMV filteredNoMV] = filterData(dataMV,rawData);

            % preprocessData normalizes the data by subtracting the mean of the 
            % metabolites for each metabolite value and dividing by the standard 
            % deviation of the metabolite. Input is the missing value data aftering 
            % filteirng (filteredMV)
            [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

            % kNN performed on the preprocessed data (scaledMV)
            [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
            % NS-kNN performed on the preprocessed data (scaledMV)
            [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

            % Return imputed data to original scale by post processing the data
            imputedData_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
            imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

            % Calculate the normalized root mean square error for each
            % imputation method
            [RMSE_kNN(iter,percentMV+1) NormalizedRMS_kNN(iter,percentMV+1)] = RMSError(imputedData_kNN_final,filteredMV,filteredNoMV);
            [RMSE_NSkNN(iter,percentMV+1) NormalizedRMS_NSkNN(iter,percentMV+1)] = RMSError(imputedData_NSkNN_final,filteredMV,filteredNoMV);
        end
        fprintf('Iteration %d: ',iter);
        toc
    end
    save(sprintf('%s_MNAR',fileName));
    clearvars -except dataSet dataFile K
end


%% MNAR-Titrate (MNAR-T) (Figure 4 and S4)

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

K = 6;

for dataSet = 1:3
    for percentMV = [10,30]
        
        % Load data
        fileName = dataFile{dataSet};
        load(fileName);
        fprintf('Working on %s, %d%% missing...\n',fileName,percentMV);
        
        for iter = 1:100
            tic
            for percentMNAR = 0:percentMV
                
                % removeDataMNART removes data using the MNAR-T modeling method
                % with inputs of the raw data (rawData), total percent
                % missingness (percentMV), and percent of all values that are
                % MNAR (percentMNAR)
                dataMV = removeDataMNART(rawData,percentMV,percentMNAR);

                % filterData function removes any samples or metabolites that have 
                % greater than 70% of the values missing. Inputs are the missing value
                % dataset (dataMV) and the original raw dataset (rawData)
                [filteredMV filteredNoMV] = filterData(dataMV,rawData);

                % preprocessData normalizes the data by subtracting the mean of the 
                % metabolites for each metabolite value and dividing by the standard 
                % deviation of the metabolite. Input is the missing value data aftering 
                % filteirng (filteredMV)
                [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

                % kNN performed on the preprocessed data (scaledMV)
                [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
                % NS-kNN performed on the preprocessed data (scaledMV)
                [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

                % Return imputed data to original scale by post processing the data
                imputedData_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
                imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

                % Calculate the normalized root mean square error for each
                % imputation method
                [RMSE_kNN(iter,percentMNAR+1) NormalizedRMS_kNN(iter,percentMNAR+1)] = RMSError(imputedData_kNN_final,filteredMV,filteredNoMV);
                [RMSE_NSkNN(iter,percentMNAR+1) NormalizedRMS_NSkNN(iter,percentMNAR+1)] = RMSError(imputedData_NSkNN_final,filteredMV,filteredNoMV);      
            end
            fprintf('Iteration %d: ',iter);
            toc
        end
        save(sprintf('%s_MNART_MV%02d',fileName,percentMV));
        clearvars -except dataSet dataFile percentMV K
    end
end

%% Mixed Missingness (MM) (Figure 5, S5-S6)

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};
plotLetterList = {'A','B','C','D'};

K = 6;

Max_percentBelowThresh_I = 70;
LowAbundThresh_II = 70;

for dataSet = 1:3
    for plotLetterNum = 1:4
        plotLetter = plotLetterList{plotLetterNum};
        
        if strcmp('A',plotLetter)
            percentMV = 10;
            percentMVlowAbund_III = 30;
        elseif strcmp('B',plotLetter)
            percentMV = 10;
            percentMVlowAbund_III = 40;
        elseif strcmp('C',plotLetter)
            percentMV = 30;
            percentMVlowAbund_III = 30;
        elseif strcmp('D',plotLetter)
            percentMV = 30;
            percentMVlowAbund_III = 40;
        end
        
        % Load data
        fileName = dataFile{dataSet};
        load(fileName);
        fprintf('Working on %s, %d%% missing, Plot %s...\n',fileName,percentMV,plotLetter);

        for iter = 1:100
            tic
            for percentBelowThresh_I = 5:Max_percentBelowThresh_I

                % removeDataMM generates missing values using the mixed
                % missingness (MM) removal method. Inputs are the raw data
                % (rawData), percent total missingness (percentMV), 
                % percentBelowThresh_I, LowAbundThresh_II, and
                % percentMVlowAbund_III. Descriptions for parameters I, II, and
                % III can be found in the Methods.
                dataMV = removeDataMM(rawData,percentMV,percentBelowThresh_I,LowAbundThresh_II,percentMVlowAbund_III);

                % filterData function removes any samples or metabolites that have 
                % greater than 70% of the values missing. Inputs are the missing value
                % dataset (dataMV) and the original raw dataset (rawData)
                [filteredMV filteredNoMV] = filterData(dataMV,rawData);

                % preprocessData normalizes the data by subtracting the mean of the 
                % metabolites for each metabolite value and dividing by the standard 
                % deviation of the metabolite. Input is the missing value data aftering 
                % filteirng (filteredMV)
                [scaledMV avgMV stddevMV] = preprocessData(filteredMV);

                % kNN performed on the preprocessed data (scaledMV)
                [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
                % NS-kNN performed on the preprocessed data (scaledMV)
                [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K);

                % Return imputed data to original scale by post processing the data
                imputedData_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
                imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

                % Calculate the normalized root mean square error for each
                % imputation method
                [RMSE_kNN(iter,percentBelowThresh_I-4) NormalizedRMS_kNN(iter,percentBelowThresh_I-4)] = RMSError(imputedData_kNN_final,filteredMV,filteredNoMV);
                [RMSE_NSkNN(iter,percentBelowThresh_I-4) NormalizedRMS_NSkNN(iter,percentBelowThresh_I-4)] = RMSError(imputedData_NSkNN_final,filteredMV,filteredNoMV);
            end
            fprintf('Iteration %d: ',iter);
            toc
        end
        save(sprintf('%s_MM_%s',fileName,plotLetter));
        clearvars -except dataSet dataFile percentMV Max_percentBelowThresh_I LowAbundThresh_II percentMVlowAbund_III K plotLetterList plotLetterNum plotLetter
    end
end

%% Sensitivity analysis (Figure S16-S19)
% ***NOTE: Takes a long time to run. Faster to split up and parallelize
% code, then recombine results.

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

K = 6;

Max_percentBelowThresh_I = 70;
Max_LowAbundThresh_II = 70;
Max_percentMVlowAbund_III = 40;

for dataSet = 1:3
    for percentMV = [10,20,30]
        
         % Load data
        fileName = dataFile{dataSet};
        load(fileName);
        fprintf('Working on %s, %d%% missing...\n',fileName,percentMV);
        
        counter = 1;
        for percentBelowThresh_I = 5:5:Max_percentBelowThresh_I
            for LowAbundThresh_II = 5:5:Max_LowAbundThresh_II
                for percentMVlowAbund_III = 5:5:Max_percentMVlowAbund_III
                    tic
                    for iter = 1:100
                        SensitivityData(counter,1) = K;
                        SensitivityData(counter,2) = percentMV;
                        SensitivityData(counter,3) = percentBelowThresh_I;
                        SensitivityData(counter,4) = LowAbundThresh_II;
                        SensitivityData(counter,5) = percentMVlowAbund_III;
                        
                        % removeDataMM generates missing values using the mixed
                        % missingness (MM) removal method. Inputs are the raw data
                        % (rawData), percent total missingness (percentMV), 
                        % percentBelowThresh_I, LowAbundThresh_II, and
                        % percentMVlowAbund_III. Descriptions for parameters I, II, and
                        % III can be found in the Methods.
                        dataMV = removeDataMM(rawData,percentMV,percentBelowThresh_I,LowAbundThresh_II,percentMVlowAbund_III);

                        % filterData function removes any samples or metabolites that have 
                        % greater than 70% of the values missing. Inputs are the missing value
                        % dataset (dataMV) and the original raw dataset (rawData)
                        [filteredMV filteredNoMV] = filterData(dataMV,rawData);

                        % preprocessData normalizes the data by subtracting the mean of the 
                        % metabolites for each metabolite value and dividing by the standard 
                        % deviation of the metabolite. Input is the missing value data aftering 
                        % filteirng (filteredMV)
                        [scaledMV avgMV stddevMV] = preprocessData(filteredMV);



                        % kNN performed on the preprocessed data (scaledMV)
                        [imputedData_kNN imputedDataWeighted_kNN] = kNNData(scaledMV,K);
                        % NS-kNN performed on the preprocessed data (scaledMV)
                        [imputedData_NSkNN imputedDataWeighted_NSkNN] = NSkNNData(scaledMV,K); % min value

                        % Return imputed data to original scale by post processing the data
                        imputedData_kNN_final = postprocessData(imputedDataWeighted_kNN, avgMV, stddevMV);
                        imputedData_NSkNN_final = postprocessData(imputedDataWeighted_NSkNN, avgMV, stddevMV);

                        % Calculate the normalized root mean square error for each
                        % imputation method
                        [RMSE_kNN NormalizedRMSE_kNN] = RMSError(imputedData_kNN_final,filteredMV,filteredNoMV);
                        [RMSE_NSkNN NormalizedRMSE_NSkNN] = RMSError(imputedData_NSkNN_final,filteredMV,filteredNoMV);


                        SensitivityData(counter,6) = NormalizedRMSE_kNN;
                        SensitivityData(counter,7) = NormalizedRMSE_NSkNN;
                        counter = counter+1;
                    end
                    fprintf('Iteration %d: ',iter);
                    toc
                end
            end
        end
        save(sprintf('%s_MM_MV%d_Sensitivity',fileName,percentMV),'SensitivityData');
    end
    clearvars -except dataSet dataFile i percentMV Max_percentBelowThresh_I Max_LowAbundThresh_II Max_percentMVlowAbund_III K
end