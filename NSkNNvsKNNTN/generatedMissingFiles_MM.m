clear;clc

dataFile = {'BacteriaData','MouseData','HumanData','UrinaryHumanData','NeuralHumanData','MicrobeData','AntibioticMouseData','TobaccoData','SMData'};

LowAbundThresh_II = 70;

for dataSet = 1:9
    load(sprintf('%s',dataFile{dataSet}));
    for percMV = [10 30]
        for percMNAR = [33 66]
            for percentMVlowAbund_III = [30 40]
                syms percentBelowThresh_I
                eqn = percMNAR == 100*(((percentBelowThresh_I)/100)*(LowAbundThresh_II/100)*(percentMVlowAbund_III/100)...
                + ((percentBelowThresh_I)/100)*(1-LowAbundThresh_II/100)*(0.5*percentMVlowAbund_III/100))/(percMV/100);
                percentBelowThresh_I = double(solve(eqn,percentBelowThresh_I));
                for i = 1:100             
                    fileName = sprintf('%s_MM_PercMV-%02d_ThreshIII-%02d_PercMNAR-%02d_rep-%03d.csv',dataFile{dataSet},percMV,percentMVlowAbund_III,percMNAR,i);

                    dataMV = removeDataMM(rawData,percMV,percentBelowThresh_I,LowAbundThresh_II,percentMVlowAbund_III);
                    dlmwrite(fileName,dataMV,'delimiter',',','precision',9);
                    clear dataMV
                end
            end
        end
    end
end
