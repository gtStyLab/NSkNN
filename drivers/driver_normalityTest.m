% Table S2: Calculation of normality in error distributions using a Kolmogorov-Smirnov test

clear;clc;

dataFile = {'BacteriaData','MouseData','HumanData'};

for dataSet = 1:3
    
    load(sprintf('%s_MCAR.mat',dataFile{dataSet}));
    
    xaxis = 0:3:percentMV;
    NormalizedRMS_kNN = NormalizedRMS_kNN(:,xaxis+1);
    NormalizedRMS_NSkNN = NormalizedRMS_NSkNN(:,xaxis+1);
    
    count_kNN = 1;
    count_NSkNN = 1;
    for i = 1:size(NormalizedRMS_kNN,2)
        if ~any(isinf(NormalizedRMS_kNN(:,i)))
            kNN_normality(count_kNN) = kstest(zscore(NormalizedRMS_kNN(:,i)));
            count_kNN = count_kNN+1;
        end
        if ~any(isinf(NormalizedRMS_NSkNN(:,i)))
            NSkNN_normality(count_NSkNN) = kstest(zscore(NormalizedRMS_NSkNN(:,i)));
            count_NSkNN = count_NSkNN+1;
        end
    end
    
    numNormal_kNN = sum(kNN_normality==0);
    numNormal_NSkNN = sum(NSkNN_normality==0);
    
    fprintf('%s Data, MCAR: %d/%d are normal using kNN\n',dataFile{dataSet}(1:end-4),numNormal_kNN,length(kNN_normality));
    fprintf('%s Data, MCAR: %d/%d are normal using NS-kNN\n',dataFile{dataSet}(1:end-4),numNormal_NSkNN,length(NSkNN_normality));
    clearvars -except dataFile dataSet
end
fprintf('--------------------\n');

for dataSet = 1:3
    
    load(sprintf('%s_MNAR.mat',dataFile{dataSet}));
    
    xaxis = 0:3:percentMV;
    NormalizedRMS_kNN = NormalizedRMS_kNN(:,xaxis+1);
    NormalizedRMS_NSkNN = NormalizedRMS_NSkNN(:,xaxis+1);
    
    count_kNN = 1;
    count_NSkNN = 1;
    for i = 1:size(NormalizedRMS_kNN,2)
        if ~any(isinf(NormalizedRMS_kNN(:,i)))
            kNN_normality(count_kNN) = kstest(zscore(NormalizedRMS_kNN(:,i)));
            count_kNN = count_kNN+1;
        end
        if ~any(isinf(NormalizedRMS_NSkNN(:,i)))
            NSkNN_normality(count_NSkNN) = kstest(zscore(NormalizedRMS_NSkNN(:,i)));
            count_NSkNN = count_NSkNN+1;
        end
    end
    
    numNormal_kNN = sum(kNN_normality==0);
    numNormal_NSkNN = sum(NSkNN_normality==0);
    
    fprintf('%s Data, MNAR: %d/%d are normal using kNN\n',dataFile{dataSet}(1:end-4),numNormal_kNN,length(kNN_normality));
    fprintf('%s Data, MNAR: %d/%d are normal using NS-kNN\n',dataFile{dataSet}(1:end-4),numNormal_NSkNN,length(NSkNN_normality));
    clearvars -except dataFile dataSet
end
fprintf('--------------------\n');

for dataSet = 1:3
    
    load(sprintf('%s_MNART_MV10.mat',dataFile{dataSet}));
    
    xaxis = 0:percentMV;
    NormalizedRMS_kNN = NormalizedRMS_kNN(:,xaxis+1);
    NormalizedRMS_NSkNN = NormalizedRMS_NSkNN(:,xaxis+1);
    
    count_kNN = 1;
    count_NSkNN = 1;
    for i = 1:size(NormalizedRMS_kNN,2)
        if ~any(isinf(NormalizedRMS_kNN(:,i)))
            kNN_normality(count_kNN) = kstest(zscore(NormalizedRMS_kNN(:,i)));
            count_kNN = count_kNN+1;
        end
        if ~any(isinf(NormalizedRMS_NSkNN(:,i)))
            NSkNN_normality(count_NSkNN) = kstest(zscore(NormalizedRMS_NSkNN(:,i)));
            count_NSkNN = count_NSkNN+1;
        end
    end
    
    numNormal_kNN = sum(kNN_normality==0);
    numNormal_NSkNN = sum(NSkNN_normality==0);
    
    fprintf('%s Data, MNAR-T, 10%% MV: %d/%d are normal using kNN\n',dataFile{dataSet}(1:end-4),numNormal_kNN,length(kNN_normality));
    fprintf('%s Data, MNAR-T, 10%% MV: %d/%d are normal using NS-kNN\n',dataFile{dataSet}(1:end-4),numNormal_NSkNN,length(NSkNN_normality));
    clearvars -except dataFile dataSet
end
fprintf('--------------------\n');

for dataSet = 1:3
    
    load(sprintf('%s_MNART_MV30.mat',dataFile{dataSet}));
    
    xaxis = 0:3:percentMV;
    NormalizedRMS_kNN = NormalizedRMS_kNN(:,xaxis+1);
    NormalizedRMS_NSkNN = NormalizedRMS_NSkNN(:,xaxis+1);
    
    count_kNN = 1;
    count_NSkNN = 1;
    for i = 1:size(NormalizedRMS_kNN,2)
        if ~any(isinf(NormalizedRMS_kNN(:,i)))
            kNN_normality(count_kNN) = kstest(zscore(NormalizedRMS_kNN(:,i)));
            count_kNN = count_kNN+1;
        end
        if ~any(isinf(NormalizedRMS_NSkNN(:,i)))
            NSkNN_normality(count_NSkNN) = kstest(zscore(NormalizedRMS_NSkNN(:,i)));
            count_NSkNN = count_NSkNN+1;
        end
    end
    
    numNormal_kNN = sum(kNN_normality==0);
    numNormal_NSkNN = sum(NSkNN_normality==0);
    
    fprintf('%s Data, MNAR-T, 30%% MV: %d/%d are normal using kNN\n',dataFile{dataSet}(1:end-4),numNormal_kNN,length(kNN_normality));
    fprintf('%s Data, MNAR-T, 30%% MV: %d/%d are normal using NS-kNN\n',dataFile{dataSet}(1:end-4),numNormal_NSkNN,length(NSkNN_normality));
    clearvars -except dataFile dataSet
end
fprintf('--------------------\n');

for dataSet = 1:3
    
    load(sprintf('%s_MM_A.mat',dataFile{dataSet}));
    
    percentMNAR = 100*(((5:percentBelowThresh_I)/100)*(LowAbundThresh_II/100)*(percentMVlowAbund_III/100)...
        + ((5:percentBelowThresh_I)/100)*(1-LowAbundThresh_II/100)*(0.5*percentMVlowAbund_III/100))/(percentMV/100);
    OneHundredPercMNAR = find(percentMNAR > 100);
    if isempty(OneHundredPercMNAR)
        OneHundredPercMNAR = length(percentMNAR)+1;
    end
    xaxis = 1:3:OneHundredPercMNAR(1)-1;
    NormalizedRMS_kNN = NormalizedRMS_kNN(:,xaxis+1);
    NormalizedRMS_NSkNN = NormalizedRMS_NSkNN(:,xaxis+1);
    
    count_kNN = 1;
    count_NSkNN = 1;
    for i = 1:size(NormalizedRMS_kNN,2)
        if ~any(isinf(NormalizedRMS_kNN(:,i)))
            kNN_normality(count_kNN) = kstest(zscore(NormalizedRMS_kNN(:,i)));
            count_kNN = count_kNN+1;
        end
        if ~any(isinf(NormalizedRMS_NSkNN(:,i)))
            NSkNN_normality(count_NSkNN) = kstest(zscore(NormalizedRMS_NSkNN(:,i)));
            count_NSkNN = count_NSkNN+1;
        end
    end
    
    numNormal_kNN = sum(kNN_normality==0);
    numNormal_NSkNN = sum(NSkNN_normality==0);
    
    fprintf('%s Data, MM, 10%% MV, III = 30%%: %d/%d are normal using kNN\n',dataFile{dataSet}(1:end-4),numNormal_kNN,length(kNN_normality));
    fprintf('%s Data, MM, 10%% MV, III = 30%%: %d/%d are normal using NS-kNN\n',dataFile{dataSet}(1:end-4),numNormal_NSkNN,length(NSkNN_normality));
    clearvars -except dataFile dataSet
end
fprintf('--------------------\n');

for dataSet = 1:3
    
    load(sprintf('%s_MM_B.mat',dataFile{dataSet}));
    
    percentMNAR = 100*(((5:percentBelowThresh_I)/100)*(LowAbundThresh_II/100)*(percentMVlowAbund_III/100)...
        + ((5:percentBelowThresh_I)/100)*(1-LowAbundThresh_II/100)*(0.5*percentMVlowAbund_III/100))/(percentMV/100);
    OneHundredPercMNAR = find(percentMNAR > 100);
    if isempty(OneHundredPercMNAR)
        OneHundredPercMNAR = length(percentMNAR)+1;
    end
    xaxis = 1:3:OneHundredPercMNAR(1)-1;
    NormalizedRMS_kNN = NormalizedRMS_kNN(:,xaxis+1);
    NormalizedRMS_NSkNN = NormalizedRMS_NSkNN(:,xaxis+1);
    
    count_kNN = 1;
    count_NSkNN = 1;
    for i = 1:size(NormalizedRMS_kNN,2)
        if ~any(isinf(NormalizedRMS_kNN(:,i)))
            kNN_normality(count_kNN) = kstest(zscore(NormalizedRMS_kNN(:,i)));
            count_kNN = count_kNN+1;
        end
        if ~any(isinf(NormalizedRMS_NSkNN(:,i)))
            NSkNN_normality(count_NSkNN) = kstest(zscore(NormalizedRMS_NSkNN(:,i)));
            count_NSkNN = count_NSkNN+1;
        end
    end
    
    numNormal_kNN = sum(kNN_normality==0);
    numNormal_NSkNN = sum(NSkNN_normality==0);
    
    fprintf('%s Data, MM, 10%% MV, III = 40%%: %d/%d are normal using kNN\n',dataFile{dataSet}(1:end-4),numNormal_kNN,length(kNN_normality));
    fprintf('%s Data, MM, 10%% MV, III = 40%%: %d/%d are normal using NS-kNN\n',dataFile{dataSet}(1:end-4),numNormal_NSkNN,length(NSkNN_normality));
    clearvars -except dataFile dataSet
end
fprintf('--------------------\n');

for dataSet = 1:3
    
    load(sprintf('%s_MM_C.mat',dataFile{dataSet}));
    
    percentMNAR = 100*(((5:percentBelowThresh_I)/100)*(LowAbundThresh_II/100)*(percentMVlowAbund_III/100)...
        + ((5:percentBelowThresh_I)/100)*(1-LowAbundThresh_II/100)*(0.5*percentMVlowAbund_III/100))/(percentMV/100);
    OneHundredPercMNAR = find(percentMNAR > 100);
    if isempty(OneHundredPercMNAR)
        OneHundredPercMNAR = length(percentMNAR)+1;
    end
    xaxis = 1:3:OneHundredPercMNAR(1)-1;
    NormalizedRMS_kNN = NormalizedRMS_kNN(:,xaxis+1);
    NormalizedRMS_NSkNN = NormalizedRMS_NSkNN(:,xaxis+1);
    
    count_kNN = 1;
    count_NSkNN = 1;
    for i = 1:size(NormalizedRMS_kNN,2)
        if ~any(isinf(NormalizedRMS_kNN(:,i)))
            kNN_normality(count_kNN) = kstest(zscore(NormalizedRMS_kNN(:,i)));
            count_kNN = count_kNN+1;
        end
        if ~any(isinf(NormalizedRMS_NSkNN(:,i)))
            NSkNN_normality(count_NSkNN) = kstest(zscore(NormalizedRMS_NSkNN(:,i)));
            count_NSkNN = count_NSkNN+1;
        end
    end
    
    numNormal_kNN = sum(kNN_normality==0);
    numNormal_NSkNN = sum(NSkNN_normality==0);
    
    fprintf('%s Data, MM, 30%% MV, III = 30%%: %d/%d are normal using kNN\n',dataFile{dataSet}(1:end-4),numNormal_kNN,length(kNN_normality));
    fprintf('%s Data, MM, 30%% MV, III = 30%%: %d/%d are normal using NS-kNN\n',dataFile{dataSet}(1:end-4),numNormal_NSkNN,length(NSkNN_normality));
    clearvars -except dataFile dataSet
end
fprintf('--------------------\n');

for dataSet = 1:3
    
    load(sprintf('%s_MM_D.mat',dataFile{dataSet}));
    
    percentMNAR = 100*(((5:percentBelowThresh_I)/100)*(LowAbundThresh_II/100)*(percentMVlowAbund_III/100)...
        + ((5:percentBelowThresh_I)/100)*(1-LowAbundThresh_II/100)*(0.5*percentMVlowAbund_III/100))/(percentMV/100);
    OneHundredPercMNAR = find(percentMNAR > 100);
    if isempty(OneHundredPercMNAR)
        OneHundredPercMNAR = length(percentMNAR)+1;
    end
    xaxis = 1:3:OneHundredPercMNAR(1)-1;
    NormalizedRMS_kNN = NormalizedRMS_kNN(:,xaxis+1);
    NormalizedRMS_NSkNN = NormalizedRMS_NSkNN(:,xaxis+1);
    
    count_kNN = 1;
    count_NSkNN = 1;
    for i = 1:size(NormalizedRMS_kNN,2)
        if ~any(isinf(NormalizedRMS_kNN(:,i)))
            kNN_normality(count_kNN) = kstest(zscore(NormalizedRMS_kNN(:,i)));
            count_kNN = count_kNN+1;
        end
        if ~any(isinf(NormalizedRMS_NSkNN(:,i)))
            NSkNN_normality(count_NSkNN) = kstest(zscore(NormalizedRMS_NSkNN(:,i)));
            count_NSkNN = count_NSkNN+1;
        end
    end
    
    numNormal_kNN = sum(kNN_normality==0);
    numNormal_NSkNN = sum(NSkNN_normality==0);
    
    fprintf('%s Data, MM, 30%% MV, III = 40%%: %d/%d are normal using kNN\n',dataFile{dataSet}(1:end-4),numNormal_kNN,length(kNN_normality));
    fprintf('%s Data, MM, 30%% MV, III = 40%%: %d/%d are normal using NS-kNN\n',dataFile{dataSet}(1:end-4),numNormal_NSkNN,length(NSkNN_normality));
    clearvars -except dataFile dataSet
end
fprintf('--------------------\n');