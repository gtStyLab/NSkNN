rm(list=ls())
cat("\014")

source("Imputation_Functions.R")
dataFile = c("BacteriaData","MouseData","HumanData")

for (dataSet in c(1,2,3)){
  dataName = dataFile[dataSet]
  for (j in c(9,15,30)){
    for (k in c(1/3,2/3)){
      mnarPer = k*j
      KnnTruncNRMSE <- vector()
      for (i in seq(1,100,length=100)){
        print(i)
        start_time <- Sys.time()
        CompleteData = read.csv(sprintf('rawData\\%s.csv',dataName),header=FALSE)
        CompleteData <- as.matrix(CompleteData)
        fileName = sprintf('MVdata\\%s\\%s_PercentMV-%02d_PercentMNAR-%02d_rep-%03d.csv',dataName,dataName,j,mnarPer,i)
        MissData = read.csv(fileName,header=FALSE)
        MissData <- as.matrix(MissData)
        MissData[is.nan(MissData)]<-NA
        
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
        
        KnnTrunc <- Nrmse(kNN_Trunc_Imp, t(MissData), t(CompleteData))
        KnnTruncNRMSE <- append(KnnTruncNRMSE,KnnTrunc)
        #names(RMSError) <- c("KNN-TN", "KNN-CR", "KNN-EU")
        rm(MissData,CompleteData,KnnTrunc)
        end_time <- Sys.time()
        print(end_time - start_time)
        print(KnnTruncNRMSE)
        
      }
      saveName = sprintf('%s_PercentMV-%02d_PercentMNAR-%02d_results.Rdata',dataName,j,mnarPer,i)
      save(list = ls(all.names = TRUE),file=saveName)
    }
  }
}
