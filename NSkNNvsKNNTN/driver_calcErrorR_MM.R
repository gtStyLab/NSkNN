rm(list=ls())
cat("\014")

source("Imputation_Functions.R")
dataFile = c("BacteriaData","MouseData","HumanData","UrinaryHumanData","NeuralHumanData","MicrobeData","AntibioticMouseData","TobaccoData","SMData")

for (dataSet in c(1:9)){
  dataName = dataFile[dataSet]
  for (percMV in c(10,30)){
    for (percMNAR in c(25,50,75)){
      for (percentMVlowAbund_III in c(30,40)){
      KnnTruncNRMSE <- vector()
        for (i in seq(1,100,length=100)){
          print(i)
          start_time <- Sys.time()
          CompleteData = read.csv(sprintf('rawData\\%s.csv',dataName),header=FALSE)
          CompleteData <- as.matrix(CompleteData)
          fileName = sprintf('MVdata\\%s_MM\\%s_MM_PercMV-%02d_ThreshIII-%02d_PercMNAR-%02d_rep-%03d.csv',dataFile[dataSet],dataFile[dataSet],percMV,percentMVlowAbund_III,percMNAR,i)
          MissData = read.csv(fileName,header=FALSE)
          MissData <- as.matrix(MissData)
          MissData[is.nan(MissData)]<-NA
          MissData <- log(MissData)
          
          # Remove rows (metabolites) with more than 75% missing data
          perc = 0.75
          NumberMissing <- apply(MissData, 2, function(x) length(which(is.na(x))))
          ind  <- which(NumberMissing >= (nrow(CompleteData)*perc))
          if(length(ind) > 0 ){
            MissData <- MissData[,-ind]
            CompleteData <- CompleteData[,-ind]
          }
          MissData <- t(MissData)
          CompleteData <- t(CompleteData)
          
          #LOD <- quantile(CompleteData, probs = 0.1)
          (kNN_Trunc_Imp <- imputeKNN(t(MissData), k=10 , distance = "truncation", perc= 0.75))
          
          kNN_Trunc_Imp <- exp(kNN_Trunc_Imp)
          MissData <- exp(MissData)
          KnnTrunc <- Nrmse(kNN_Trunc_Imp, t(MissData), t(CompleteData))
          KnnTruncNRMSE <- append(KnnTruncNRMSE,KnnTrunc)
          #names(RMSError) <- c("KNN-TN", "KNN-CR", "KNN-EU")
          rm(MissData,CompleteData,KnnTrunc)
          end_time <- Sys.time()
          print(end_time - start_time)
          print(KnnTruncNRMSE)
          
        }
        saveName = sprintf('%s_MM_PercMV-%02d_ThreshIII-%02d_PercMNAR-%02d_results.Rdata',dataName,percMV,percentMVlowAbund_III,percMNAR)
        save(list = ls(all.names = TRUE),file=saveName)
      }
    }
  }
}
