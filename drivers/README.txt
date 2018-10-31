driver_GeneratePlotData.m generates the data for all figures and tables in the manuscript and supplementary, except for Figure 6, S8-S15 (see folder NSkNNvsKNNTN). Data for Figure S7 data is also not generated, as it just uses the raw data.

driver_normalityTest.m tests the normality of the error distributions of some of the data generated from driver_GeneratePlotData.m using as Kolmogorov-Smirnov test.

driver_NSkNN.m is an example of how NSkNN would be run on a metabolomics dataset with missing values, stored as a CSV file.

driver_PlotPaperFigures.m plots all the figures in the manuscript and supplementary using the data generated from driver_GeneratePlotData.m and NSkNNvsKNNTN. Legacy result files for recreating exact figures in the paper are found in the results folder.