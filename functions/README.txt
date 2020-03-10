filterData.m filters the samples and metabolites that have >70% of values missing (by default).

kNNData.m imputes data using the original kNN algorithm with Euclidean distance as a distance metric that skips nearest neighbors with missing values in the same metabolite location as the value in the target sample being imputed.

kNNData_PearsonCorr.m  imputes data using the original kNN algorithm with Pearson correlation as a distance metric that skips nearest neighbors with missing values in the same metabolite location as the value in the target sample being imputed.

NSkNNData.m imputes data using No Skip-kNN that replaces missing values in the imputation calculation with the minimum abundance of the metabolite.

NSkNNData_HM.m imputes data using No Skip-kNN that replaces missing values in the imputation calculation with the half minimum abundance of the metabolite.

NSkNNData_Zero.m imputes data using No Skip-kNN that replaces missing values in the imputation calculation with zero.

postprocessData.m reverts the imputed data (if autoscaled) back to its original scale.

preprocessData.m autoscales the data by subtracting the mean and dividing by the standard deviation of each metabolite.

removeDataMCAR.m generates missing values for a complete dataset by the MCAR modeling method.

removeDataMM.m generates missing values for a complete dataset by the MM modeling method.

removeDataMNAR.m generates missing values for a complete dataset by the MNAR modeling method.

removedDataMNAR_LODparam.m generates missing values for a complete dataset by the MNAR modeling method, with an LOD percentage threshold as an input. This is used for Figure S1 in the Supplementary.

removeDataMNART.m generates missing values for a complete dataset by the MNAR-T modeling method.

RMSError.m calculates the root mean squared error and normalized root mean squared error of the imputed data.