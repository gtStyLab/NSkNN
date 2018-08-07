%% MNAR testing LOD parameter (Figure S1)

clear;clc;
dataFile = {'BacteriaData','MouseData','HumanData'};
figureLetter = {'a.','b.','c.'};
annotPosition = [0.1 0.98 0 0; 0.54 0.98 0 0; 0.1 0.51 0 0];
figure('units','normalized','outerposition',[0 0 0.75 1]);

for dataSet = 1:3
    for LOD = [20 40 60 80]
        subplot(2,2,dataSet);
        figAnnot = annotation('textbox', annotPosition(dataSet,:), 'string', figureLetter{dataSet});
        figAnnot.FontSize = 18;
        figAnnot.FontWeight = 'bold';

        load(sprintf('%s_MNAR_LOD%d.mat',dataFile{dataSet},LOD));

        hold on;

        meanNRMSE_kNN = mean(NormalizedRMS_kNN);
        stdNRMSE_kNN = std(NormalizedRMS_kNN);
        meanNRMSE_NSkNN = mean(NormalizedRMS_NSkNN);
        stdNRMSE_NSkNN = std(NormalizedRMS_NSkNN);
        
        % Plot all points or a select amount of points
        %xaxis = 0:percentMV;
        xaxis = 0:3:percentMV;

        a = errorbar(xaxis,meanNRMSE_kNN(xaxis+1),stdNRMSE_kNN(xaxis+1));
        b = errorbar(xaxis,meanNRMSE_NSkNN(xaxis+1),stdNRMSE_NSkNN(xaxis+1));

        a.Color = 'b';
        b.Color = 'g';

        if LOD == 20
            a.Marker = 'p';
            a.LineWidth = 1.5;
            b.Marker = 'p';
            b.LineWidth = 1.5;
            b.LineStyle = '--';
        elseif LOD == 40
            a.Marker = 's';
            a.LineWidth = 1.5;
            b.Marker = 's';
            b.LineWidth = 1.5;
            b.LineStyle = '--';
        elseif LOD == 60
            a.Marker = '*';
            a.LineWidth = 1.5;
            b.Marker = '*';
            b.LineWidth = 1.5;
            b.LineStyle = '--';
        elseif LOD == 80
            a.Marker = 'o';
            a.LineWidth = 1.5;
            b.Marker = 'o';
            b.LineWidth = 1.5;
            b.LineStyle = '--';
        end
        
        a.MarkerSize = 11;
        b.MarkerSize = 11;

        subplot(2,2,dataSet);
        figAnnot = annotation('textbox', annotPosition(dataSet,:), 'string', figureLetter{dataSet});
        figAnnot.FontSize = 18;
        figAnnot.FontWeight = 'bold';

        set(gca,'fontsize',14);
        title(sprintf('MNAR %s Dataset',dataFile{dataSet}(1:end-4)),'FontSize',18);
        xlabel('% Total Missing');
        ylabel('Normalized Root Mean Square Error');
        ylim([-inf 0.65]);
        if dataSet == 1  
            lgd1 = line(nan,nan,'Linestyle','-','LineWidth',1,'Marker','p','Color','b');
            lgd2 = line(nan,nan,'Linestyle','--','Marker','p','Color','g');
            lgd3 = line(nan,nan,'Linestyle','-','LineWidth',1,'Marker','s','Color','b');
            lgd4 = line(nan,nan,'Linestyle','--','Marker','s','Color','g');
            lgd5 = line(nan,nan,'Linestyle','-','LineWidth',1,'Marker','*','Color','b');
            lgd6 = line(nan,nan,'Linestyle','--','Marker','*','Color','g');
            lgd7 = line(nan,nan,'Linestyle','-','LineWidth',1,'Marker','o','Color','b');
            lgd8 = line(nan,nan,'Linestyle','--','Marker','o','Color','g');
            lgd = legend([lgd1,lgd2,lgd3,lgd4,lgd5,lgd6,lgd7,lgd8],{'kNN - 20% LOD Threshold','NS-kNN - 20% LOD Threshold','kNN - 40% LOD Threshold','NS-kNN - 40% LOD Threshold','kNN - 60% LOD Threshold','NS-kNN - 60% LOD Threshold','kNN - 80% LOD Threshold','NS-kNN - 80% LOD Threshold'},'Location','northwest');
            lgd.FontSize = 10;
        end
    end
end


%% MNAR-Titrate (MNAR-T) comparison of kNN, NSkNN, NSkNN_HM, NSkNN_Zero (Figure S2)

clear;clc;
dataFile = {'BacteriaData','MouseData','HumanData'};
figureLetter = {'a.','b.','c.'};
annotPosition = [0.1 0.98 0 0; 0.54 0.98 0 0; 0.1 0.51 0 0];

for percentMV = [30]
    figure('units','normalized','outerposition',[0 0 0.75 1]);
    for dataSet = [1:3]

        subplot(2,2,dataSet);
        figAnnot = annotation('textbox', annotPosition(dataSet,:), 'string', figureLetter{dataSet});
        figAnnot.FontSize = 18;
        figAnnot.FontWeight = 'bold';

        load(sprintf('%s_MNART_MV%02d_compareNSkNNReplacement.mat',dataFile{dataSet},percentMV));

        hold on;

        meanNRMSE_kNN = mean(NormalizedRMS_kNN);
        stdNRMSE_kNN = std(NormalizedRMS_kNN);
        meanNRMSE_NSkNN = mean(NormalizedRMS_NSkNN);
        stdNRMSE_NSkNN = std(NormalizedRMS_NSkNN);
        meanNRMSE_NSkNN_HM = mean(NormalizedRMS_NSkNN_HM);
        stdNRMSE_NSkNN_HM = std(NormalizedRMS_NSkNN_HM);
        meanNRMSE_NSkNN_zero = mean(NormalizedRMS_NSkNN_zero);
        stdNRMSE_NSkNN_zero = std(NormalizedRMS_NSkNN_zero);
        
        % Plot all points or a select amount of points
        %xaxis = 0:percentMNAR;
        xaxis = 0:3:percentMNAR;
        
        a = errorbar(xaxis,meanNRMSE_kNN(xaxis+1),stdNRMSE_kNN(xaxis+1));
        b = errorbar(xaxis,meanNRMSE_NSkNN(xaxis+1),stdNRMSE_NSkNN(xaxis+1));
        c = errorbar(xaxis,meanNRMSE_NSkNN_HM(xaxis+1),stdNRMSE_NSkNN_HM(xaxis+1));
        d = errorbar(xaxis,meanNRMSE_NSkNN_zero(xaxis+1),stdNRMSE_NSkNN_zero(xaxis+1));
        
        a.Color = 'b';
        a.LineStyle = '-';
        a.LineWidth = 1.5;
        b.Color = 'g';
        b.LineStyle = '--';
        b.LineWidth = 1.5;
        c.Color = 'r';
        c.LineStyle = ':';
        c.LineWidth = 1.5;
        d.Color = 'k';
        d.LineStyle = '-.';
        d.LineWidth = 1.5;

        set(gca,'fontsize',14);
        title(sprintf('MNAR-T %s Dataset',dataFile{dataSet}(1:end-4)),'FontSize',18);
        xlabel('% missing values that are MNAR');
        ylabel('Normalized Root Mean Square Error');
        set(gca,'Xtick',[0:(max(percentMNAR)/5):max(percentMNAR)],'XTickLabel',[1:6])
        set(gca,'XtickLabel',{'0' '20' '40' '60' '80' '100'})
        if dataSet == 1
            lgd = legend('kNN','NS-kNN','NS-kNN HM','NS-kNN Zero','Location','southwest');
            lgd.FontSize = 14;
        end
    end
end


%% MCAR Best K value (Figure S3)

clear;clc;
dataFile = {'BacteriaData','MouseData','HumanData'};
figureLetter = {'a.','b.','c.'};
annotPosition = [0.1 0.98 0 0; 0.54 0.98 0 0; 0.1 0.51 0 0];
figure('units','normalized','outerposition',[0 0 0.75 1]);

for dataSet = 1:3
    subplot(2,2,dataSet);
    figAnnot = annotation('textbox', annotPosition(dataSet,:), 'string', figureLetter{dataSet});
    figAnnot.FontSize = 18;
    figAnnot.FontWeight = 'bold';
    
    load(sprintf('%s_MCAR_K.mat',dataFile{dataSet}));
    
    hold on;
    
    xaxis = 1:MaxK;
        
    % Plot line graphs
    a = errorbar(xaxis,mean(NormalizedRMS_kNN),std(NormalizedRMS_kNN));
    b = errorbar(xaxis,mean(NormalizedRMS_NSkNN),std(NormalizedRMS_NSkNN));
    
    % Plot significance
    max_yaxis = max([mean(NormalizedRMS_kNN) + std(NormalizedRMS_kNN); mean(NormalizedRMS_NSkNN) + std(NormalizedRMS_NSkNN)]) + max([mean(NormalizedRMS_kNN) + std(NormalizedRMS_kNN); mean(NormalizedRMS_NSkNN) + std(NormalizedRMS_NSkNN)])/40;
    max_yaxis = max_yaxis(xaxis);
    significant = ttest2(NormalizedRMS_kNN,NormalizedRMS_NSkNN);
    significant(isnan(significant)) = 0;
    plot(xaxis(logical(significant(xaxis))),max_yaxis(logical(significant(xaxis))),'k*');
    
    a.Color = 'b';
    a.LineStyle = '-';
    a.LineWidth = 1.5;
    b.Color = 'g';
    b.LineStyle = '--';
    b.LineWidth = 1.5;

    set(gca,'fontsize',14);
    title(sprintf('MCAR %s Dataset',dataFile{dataSet}(1:end-4)));
    xlabel('# of Nearest Neighbors (k)');
    ylabel('Normalized Root Mean Square Error');
    if dataSet == 1
        lgd = legend('kNN','NS-kNN','Location','northwest');
        lgd.FontSize = 14;
    end
    xlim([1 K]);
end

%% MCAR increasing % MV (Figure 2)

clear;clc;
dataFile = {'BacteriaData','MouseData','HumanData'};
figureLetter = {'a.','b.','c.'};
annotPosition = [0.1 0.98 0 0; 0.54 0.98 0 0; 0.1 0.51 0 0];
figure('units','normalized','outerposition',[0 0 0.26 0.43]);

for dataSet = 1:3
    subplot(2,2,dataSet);
    figAnnot = annotation('textbox', annotPosition(dataSet,:), 'string', figureLetter{dataSet});
    figAnnot.FontSize = 11;
    figAnnot.FontWeight = 'bold';
    
    load(sprintf('%s_MCAR.mat',dataFile{dataSet}));

    hold on;

    meanNRMSE_kNN = mean(NormalizedRMS_kNN);
    stdNRMSE_kNN = std(NormalizedRMS_kNN);
    meanNRMSE_NSkNN = mean(NormalizedRMS_NSkNN);
    stdNRMSE_NSkNN = std(NormalizedRMS_NSkNN);
    
    % Plot all points or a select amount of points
    %xaxis = 0:percentMV;
    xaxis = 0:3:percentMV;
        
    % Plot line graphs
    a = errorbar(xaxis,meanNRMSE_kNN(xaxis+1),stdNRMSE_kNN(xaxis+1));
    b = errorbar(xaxis,meanNRMSE_NSkNN(xaxis+1),stdNRMSE_NSkNN(xaxis+1));
    
    % Plot significance
    max_yaxis = max([mean(NormalizedRMS_kNN) + std(NormalizedRMS_kNN); mean(NormalizedRMS_NSkNN) + std(NormalizedRMS_NSkNN)]) + max([mean(NormalizedRMS_kNN) + std(NormalizedRMS_kNN); mean(NormalizedRMS_NSkNN) + std(NormalizedRMS_NSkNN)])/40;
    max_yaxis = max_yaxis(xaxis+1);
    significant = ttest2(NormalizedRMS_kNN,NormalizedRMS_NSkNN);
    significant(isnan(significant)) = 0;
    plot(xaxis(logical(significant(xaxis+1))),max_yaxis(logical(significant(xaxis+1))),'k*');
    
    a.Color = 'b';
    a.LineStyle = '-';
    a.LineWidth = 1;
    b.Color = 'g';
    b.LineStyle = '--';
    b.LineWidth = 1;
    
    set(gca,'fontsize',9);
    title(sprintf('MCAR %s Dataset',dataFile{dataSet}(1:end-4)),'FontSize',10);
    xlabel('% Total Missing');
    ylabel('Normalized Root Mean Square Error');
    if dataSet == 1
        lgd = legend('kNN','NS-kNN','Location','northwest');
        lgd.FontSize = 9;
    end
end

%% MNAR (Figure 3)

clear;clc;
dataFile = {'BacteriaData','MouseData','HumanData'};
figureLetter = {'a.','b.','c.'};
annotPosition = [0.1 0.98 0 0; 0.54 0.98 0 0; 0.1 0.51 0 0];
figure('units','normalized','outerposition',[0 0 0.26 0.43]);

for dataSet = 1:3
    
    subplot(2,2,dataSet);
    figAnnot = annotation('textbox', annotPosition(dataSet,:), 'string', figureLetter{dataSet});
    figAnnot.FontSize = 11;
    figAnnot.FontWeight = 'bold';
    
    load(sprintf('%s_MNAR.mat',dataFile{dataSet}));

    hold on;
    
    meanNRMSE_kNN = mean(NormalizedRMS_kNN);
    stdNRMSE_kNN = std(NormalizedRMS_kNN);
    meanNRMSE_NSkNN = mean(NormalizedRMS_NSkNN);
    stdNRMSE_NSkNN = std(NormalizedRMS_NSkNN);
    
    % Plot all points or a select amount of points
    %xaxis = 0:percentMV;
    xaxis = 0:3:percentMV;
    
    a = errorbar(xaxis,meanNRMSE_kNN(xaxis+1),stdNRMSE_kNN(xaxis+1));
    b = errorbar(xaxis,meanNRMSE_NSkNN(xaxis+1),stdNRMSE_NSkNN(xaxis+1));
    
    % Plot significance
    max_yaxis = max([mean(NormalizedRMS_kNN) + std(NormalizedRMS_kNN); mean(NormalizedRMS_NSkNN) + std(NormalizedRMS_NSkNN)]) + max([mean(NormalizedRMS_kNN) + std(NormalizedRMS_kNN); mean(NormalizedRMS_NSkNN) + std(NormalizedRMS_NSkNN)])/40;
    max_yaxis = max_yaxis(xaxis+1);
    significant = ttest2(NormalizedRMS_kNN,NormalizedRMS_NSkNN);
    significant(isnan(significant)) = 0;
    plot(xaxis(logical(significant(xaxis+1))),max_yaxis(logical(significant(xaxis+1))),'k*');
    
    a.Color = 'b';
    a.LineStyle = '-';
    a.LineWidth = 1;
    b.Color = 'g';
    b.LineStyle = '--';
    b.LineWidth = 1;

    set(gca,'fontsize',9);
    title(sprintf('MNAR %s Dataset',dataFile{dataSet}(1:end-4)),'FontSize',10);
    xlabel('% Total Missing');
    ylabel('Normalized Root Mean Square Error');
    if dataSet == 1
        lgd = legend('kNN','NS-kNN','Location','northwest');
        lgd.FontSize = 9;
    end
end

%% MNAR-Titrate (MNAR-T) (Figure 4 and Figure S4)

clear;clc;
dataFile = {'BacteriaData','MouseData','HumanData'};
figureLetter = {'a.','b.','c.'};
annotPosition = [0.1 0.98 0 0; 0.54 0.98 0 0; 0.1 0.51 0 0];

for percentMV = [10, 30]
    figure('units','normalized','outerposition',[0 0 0.26 0.43]);
    for dataSet = 1:3

        subplot(2,2,dataSet);
        figAnnot = annotation('textbox', annotPosition(dataSet,:), 'string', figureLetter{dataSet});
        figAnnot.FontSize = 11;
        figAnnot.FontWeight = 'bold';

        load(sprintf('%s_MNART_MV%02d.mat',dataFile{dataSet},percentMV));

        hold on;

        meanNRMSE_kNN = mean(NormalizedRMS_kNN);
        stdNRMSE_kNN = std(NormalizedRMS_kNN);
        meanNRMSE_NSkNN = mean(NormalizedRMS_NSkNN);
        stdNRMSE_NSkNN = std(NormalizedRMS_NSkNN);

        % Percentage at which kNN and NS-kNN error plots intersect. Uses InterX
        % function by NS 2010 found on MathWorks File Exchange.
        %{
        intersection = InterX([0:percentMNAR;meanNRMSE_kNN],[0:percentMNAR;meanNRMSE_NSkNN]);
        if isempty(intersection)
            0 % If no intersection
        else
            intersection(1,end)/percentMV
        end
        %}
        
        % Plot all points or a select amount of points
        %xaxis = 0:percentMNAR; % Used for Figure S4
        xaxis = 0:3:percentMNAR; % Used for Figure 4
        
        a = errorbar(xaxis,meanNRMSE_kNN(xaxis+1),stdNRMSE_kNN(xaxis+1));
        b = errorbar(xaxis,meanNRMSE_NSkNN(xaxis+1),stdNRMSE_NSkNN(xaxis+1));
        
        % Plot significance       
        max_yaxis = max([mean(NormalizedRMS_kNN) + std(NormalizedRMS_kNN); mean(NormalizedRMS_NSkNN) + std(NormalizedRMS_NSkNN)]) + max([mean(NormalizedRMS_kNN) + std(NormalizedRMS_kNN); mean(NormalizedRMS_NSkNN) + std(NormalizedRMS_NSkNN)])/40;
        max_yaxis = max_yaxis(xaxis+1);
        significant = ttest2(NormalizedRMS_kNN,NormalizedRMS_NSkNN);
        significant(isnan(significant)) = 0;
        plot(xaxis(logical(significant(xaxis+1))),max_yaxis(logical(significant(xaxis+1))),'k*');
        
        a.Color = 'b';
        a.LineStyle = '-';
        a.LineWidth = 1;
        b.Color = 'g';
        b.LineStyle = '--';
        b.LineWidth = 1;

        set(gca,'fontsize',9);
        title(sprintf('MNAR-T %s Dataset',dataFile{dataSet}(1:end-4)),'FontSize',10);
        xlabel('% missing values that are MNAR');
        ylabel('Normalized Root Mean Square Error');
        set(gca,'Xtick',[0:(max(percentMNAR)/5):max(percentMNAR)],'XTickLabel',[1:6])
        set(gca,'XtickLabel',{'0' '20' '40' '60' '80' '100'})
        if dataSet == 1
            lgd = legend('kNN','NS-kNN','Location','northwest');
            lgd.FontSize = 9;
        end
    end
end

%% Mixed Missingness (MM) (Figure 5, Figure S5, and Figure S6) 

clear;clc;
dataFile = {'BacteriaData','MouseData','HumanData'};
figureLetter = {'a.','b.','c.','d.'};
annotPosition = [0.1 0.98 0 0; 0.54 0.98 0 0; 0.1 0.51 0 0; 0.54 0.51 0 0];

for dataSet = 1:3
    figure('units','normalized','outerposition',[0 0 0.26 0.43]);
    
    for figureNum = 1:3
        load(sprintf('%s_MM_%s.mat',dataFile{dataSet},figureLetter{figureNum}(1:end-1)));
        subplot(2,2,figureNum);
        figAnnot = annotation('textbox', annotPosition(figureNum,:), 'string', figureLetter{figureNum});
        figAnnot.FontSize = 11;
        figAnnot.FontWeight = 'bold';
        hold on;

        percentMNAR = 100*(((5:percentBelowThresh_I)/100)*(LowAbundThresh_II/100)*(percentMVlowAbund_III/100)...
            + ((5:percentBelowThresh_I)/100)*(1-LowAbundThresh_II/100)*(0.5*percentMVlowAbund_III/100))/(percentMV/100);
        OneHundredPercMNAR = find(percentMNAR > 100);
        if isempty(OneHundredPercMNAR)
            OneHundredPercMNAR = length(percentMNAR)+1;
        end

        mean_NormalizedRMS_kNN = mean(NormalizedRMS_kNN);
        std_NormalizedRMS_kNN = std(NormalizedRMS_kNN);
        mean_NormalizedRMS_NSkNN = mean(NormalizedRMS_NSkNN);
        std_NormalizedRMS_NSkNN = std(NormalizedRMS_NSkNN);
        
        % Percentage at which kNN and NS-kNN error plots intersect. Uses InterX
        % function by NS 2010 found on MathWorks File Exchange.
        intersection = InterX([percentMNAR(1:OneHundredPercMNAR(1)-1);mean_NormalizedRMS_kNN(1:OneHundredPercMNAR(1)-1)],[percentMNAR(1:OneHundredPercMNAR(1)-1);mean_NormalizedRMS_NSkNN(1:OneHundredPercMNAR(1)-1)]);
        if isempty(intersection)
            0 % If no intersection
        else
            intersection(1,end)
        end
        
        % Plot all points or a select amount of points
        %xaxis = 1:OneHundredPercMNAR(1)-1;
        xaxis = 1:3:OneHundredPercMNAR(1)-1;
        
        a = errorbar(percentMNAR(xaxis),mean_NormalizedRMS_kNN(xaxis),std_NormalizedRMS_kNN(xaxis));
        b = errorbar(percentMNAR(xaxis),mean_NormalizedRMS_NSkNN(xaxis),std_NormalizedRMS_NSkNN(xaxis));
        
        % Plot significance       
        max_yaxis = max([mean_NormalizedRMS_kNN + std_NormalizedRMS_kNN; mean_NormalizedRMS_NSkNN + std_NormalizedRMS_NSkNN]) + max([mean_NormalizedRMS_kNN + std_NormalizedRMS_kNN; mean_NormalizedRMS_NSkNN + std_NormalizedRMS_NSkNN])/40;
        max_yaxis = max_yaxis(xaxis);
        significant = ttest2(NormalizedRMS_kNN,NormalizedRMS_NSkNN);
        significant(isnan(significant)) = 0;
        plot(percentMNAR(xaxis(logical(significant(xaxis)))),max_yaxis(logical(significant(xaxis))),'k*');
        
        a.Color = 'b';
        a.LineStyle = '-';
        a.LineWidth = 1;
        b.Color = 'g';
        b.LineStyle = '--';
        b.LineWidth = 1;
        
        set(gca,'fontsize',9);
        title(sprintf('Percent Missing = %d%%, III = %d%%',percentMV,percentMVlowAbund_III),'FontSize',10);
        xlabel('% missing values that are MNAR');
        ylabel('Normalized Root Mean Square Error');
        if figureNum == 1
            lgd = legend('kNN','NS-kNN','Location','northwest');
            lgd.FontSize = 9;
        end
    end
end

%% NSkNN vs. kNN-TN (Figure 6, S7, S8)

clear;clc;
dataFile = {'BacteriaData','MouseData','HumanData'};

figureLetter = {'a.','b.','c.','d.','e.','f.'};
annotPosition = [0.1 0.97 0 0; 0.39 0.97 0 0; 0.67 0.97 0 0;...
    0.1 0.51 0 0; 0.39 0.51 0 0; 0.67 0.51 0 0];

for dataSet = 1:3
    figure('units','normalized','outerposition',[0 0 0.26 0.43]);
    allData_matlab = [];
    allData_csv = [];
    allData = [];
    subplot_count = 1;
    for percMNAR = [1/3 2/3]
        for percMV = [9 15 30]
            hold on;
            fileName = sprintf('%s_PercentMV-%02d_PercentMNAR-%02d_results.mat',dataFile{dataSet},percMV,percMV*percMNAR);
            csvName = sprintf('%s_PercentMV-%02d_PercentMNAR-%02d_results.csv',dataFile{dataSet},percMV,percMV*percMNAR);
            matlabData = load(fileName);
            csvData = csvread(csvName,1,1);

            allData = [matlabData.NormalizedRMS_NSkNN_compiled csvData];
            
            labelnames_allData = {'NS-kNN','kNN-TN'};
            subplot(2,3,subplot_count);
            figAnnot = annotation('textbox', annotPosition(subplot_count,:), 'string', figureLetter{subplot_count});
            figAnnot.FontSize = 9;
            figAnnot.FontWeight = 'bold';
            a = boxplot(allData,labelnames_allData);
            set(a,{'LineWidth'},{1});
                        
            set(gca,'fontsize',9);
            ylabel('Normalized RMSE');
            title({sprintf('%s Data',dataFile{dataSet}(1:end-4)),sprintf('%d%% MV, ^{%d}/_{%d} MNAR',percMV,percMNAR*3,3)},'FontSize',10);
            subplot_count = subplot_count + 1;
            ylim([0 2]);
            
            % Plot significance
            significance = ttest(allData(:,1),allData(:,2));
            if logical(significance)
                xt = get(gca, 'XTick');
                hold on;
                plot(xt([1 2]), [1 1]*min(min(allData))*0.8, '-k',  mean(xt([1 2])), min(min(allData))*0.55, '*k')
            end
        end
    end
end

%% Plot sensitivity (Figure S9, S10, S11)

clear; clc;
dataFile = {'BacteriaData','MouseData','HumanData'};
figureLetter = {'a.','b.','c.','d.','e.','f.'};
annotPosition = [0.1 0.98 0 0; 0.52 0.98 0 0;...
    0.1 0.66 0 0; 0.52 0.68 0 0;...
    0.1 0.36 0 0; 0.52 0.38 0 0];

for dataSet = 1:3
    figure('units','normalized','outerposition',[0 0 0.62 1]);
    hold on;
    
    subplot_count = 1;
    for percMV = [10 20 30]
        load(sprintf('%s_MM_MV%d_Sensitivity',dataFile{dataSet},percMV));

        for i = 1:2156
            dataStart = (i-1)*100+1;
            ErrorAvg(i,:) = mean(SensitivityData(dataStart:dataStart+99,:));
            ErrorStd(i,:) = std(SensitivityData(dataStart:dataStart+99,6:7));
            ErrorAvgStd = [ErrorAvg ErrorStd];
        end
        
        load mycmap;
        percentBelowThresh_I = ErrorAvgStd(:,3);
        LowAbundThresh_II = ErrorAvgStd(:,4);
        percentMVlowAbund_III = ErrorAvgStd(:,5);
        ErrorDiff = ErrorAvgStd(:,7)-ErrorAvgStd(:,6);

        subplot(3,2,subplot_count);
        figAnnot = annotation('textbox', annotPosition(subplot_count,:), 'string', figureLetter{subplot_count});
        figAnnot.FontSize = 18;
        figAnnot.FontWeight = 'bold';
            
        scatter3(percentBelowThresh_I,LowAbundThresh_II,percentMVlowAbund_III,20,ErrorDiff,'filled')
        xlabel('% Met. w/ MNAR (I)')
        ylabel('% Met. considered Low Abund. (II)')
        zlabel('% MNAR in Low Abund. (III)');
        xlim([0 100]);
        ylim([0 100]);
        %title(sprintf('kNN vs. NS-kNN - %s Data',dataFile{dataSet}(1:end-4)));

        cb = colorbar;
        cb.Label.String = 'NS-kNN NRMSE - kNN NRMSE';
        ax = gca;
        colormap(ax,mycmap)
        caxis([-max(abs(ErrorDiff)) max(abs(ErrorDiff))]);
        subplot_count = subplot_count + 1;
        
        % -------------
        

        percentBelowThresh_I = ErrorAvgStd(:,3);
        LowAbundThresh_II = ErrorAvgStd(:,4);
        percentMVlowAbund_III = ErrorAvgStd(:,5);
        percentMV = ErrorAvg(1,2);

        DiffErrorAvg = ErrorAvgStd(:,7)-ErrorAvgStd(:,6);
        DiffErrorStd = sqrt(ErrorAvgStd(:,9).^2+ErrorAvgStd(:,8).^2);
        
        
        percentMNAR = 100*((percentBelowThresh_I/100).*(LowAbundThresh_II/100).*(percentMVlowAbund_III/100)...
        + (percentBelowThresh_I/100).*(1-LowAbundThresh_II/100).*(0.5*percentMVlowAbund_III/100))./(percentMV/100);

        sorteverything = [percentMNAR DiffErrorAvg DiffErrorStd];
        sorteverything = sortrows(sorteverything);
        
        fewpoints = [];
        for percentThresh = [0 10 20 30 40 50 60 70 80 90 100]
            tempLoc = find(sorteverything(:,1) > percentThresh);
            if ~isempty(tempLoc)
                if tempLoc(1) == 1
                    fewpoints = [fewpoints tempLoc(1)];  
                else
                    fewpoints = [fewpoints tempLoc(1)-1];  
                end
            end
        end
        
        sortfewpoints = sorteverything(fewpoints,:);
        
        
        subplot(3,2,subplot_count);
        figAnnot = annotation('textbox', annotPosition(subplot_count,:), 'string', figureLetter{subplot_count});
        figAnnot.FontSize = 18;
        figAnnot.FontWeight = 'bold';
            
        hold on;
        a = errorbar(sortfewpoints(:,1),sortfewpoints(:,2),sortfewpoints(:,3));
        %a = errorbar(sorteverything(:,1),sorteverything(:,2),sorteverything(:,3));
        b = plot(0:1:100,zeros(1,101));
        xlim([min(sortfewpoints(:,1)) max(sortfewpoints(:,1))]);
        
        a.Color = 'b';
        a.LineWidth = 1.5;
        b.Color = 'r';
        b.LineWidth = 1.5;
        
        %title(sprintf('kNN vs. NS-kNN - %s Data',dataFile{dataSet}(1:end-4)));
        xlabel('% missing values that are MNAR');
        ylabel('NS-kNN NRMSE - kNN NRMSE');
        subplot_count = subplot_count + 1;
        
        % Percentage at which kNN and NS-kNN error plots intersect. Uses InterX
        % function by NS 2010 found on MathWorks File Exchange.
        %{
        intersection = InterX([sorteverything(:,1)';sorteverything(:,2)'],[(0:1:100);zeros(1,101)]);
        if isempty(intersection)
            0 % If no intesrection
        else
            intersection(1,end)
        end
        %}
    end
end

%% Plot sensitivity (Figure S12)

clear; clc;
figureLetter = {'a.','b.'};
annotPosition = [0.06 0.98 0 0; 0.49 0.98 0 0];

load('MouseData_MM_MV30_Sensitivity');

for i = 1:2156
    dataStart = (i-1)*100+1;
    ErrorAvg(i,:) = mean(SensitivityData(dataStart:dataStart+99,:));
    ErrorStd(i,:) = std(SensitivityData(dataStart:dataStart+99,6:7));
    ErrorAvgStd = [ErrorAvg ErrorStd];
end
        
percentBelowThresh_I = ErrorAvgStd(:,3);
LowAbundThresh_II = ErrorAvgStd(:,4);
percentMVlowAbund_III = ErrorAvgStd(:,5);
percentMV = ErrorAvg(1,2);

DiffErrorAvg = ErrorAvgStd(:,7)-ErrorAvgStd(:,6);
DiffErrorStd = sqrt(ErrorAvgStd(:,9).^2+ErrorAvgStd(:,8).^2);


percentMNAR = 100*((percentBelowThresh_I/100).*(LowAbundThresh_II/100).*(percentMVlowAbund_III/100)...
+ (percentBelowThresh_I/100).*(1-LowAbundThresh_II/100).*(0.5*percentMVlowAbund_III/100))./(percentMV/100);

sorteverything = [percentMNAR DiffErrorAvg DiffErrorStd];
sorteverything = sortrows(sorteverything);

fewpoints = [];
for percentThresh = [0 10 20 30 40 50 60 70 80 90 100]
    tempLoc = find(sorteverything(:,1) > percentThresh);
    if ~isempty(tempLoc)
        if tempLoc(1) == 1
            fewpoints = [fewpoints tempLoc(1)];  
        else
            fewpoints = [fewpoints tempLoc(1)-1];  
        end
    end
end

sortfewpoints = sorteverything(fewpoints,:);

figure('units','normalized','outerposition',[0 0 0.75 0.75]);

subplot(1,2,1);
hold on;
figAnnot = annotation('textbox', annotPosition(1,:), 'string', figureLetter{1});
figAnnot.FontSize = 18;
figAnnot.FontWeight = 'bold';
a = errorbar(sorteverything(:,1),sorteverything(:,2),sorteverything(:,3));
%plot(sorteverything(:,1),sorteverything(:,2),'m','LineWidth',2);
b = plot(0:1:100,zeros(1,101));
xlim([min(sortfewpoints(:,1)) max(sortfewpoints(:,1))]);
set(gca,'fontsize',14);
xlabel('% missing values that are MNAR');
ylabel('NS-kNN NRMSE - kNN NRMSE');

subplot(1,2,2);
hold on;
figAnnot = annotation('textbox', annotPosition(2,:), 'string', figureLetter{2});
figAnnot.FontSize = 18;
figAnnot.FontWeight = 'bold';
c = errorbar(sortfewpoints(:,1),sortfewpoints(:,2),sortfewpoints(:,3));
d = plot(0:1:100,zeros(1,101));
xlim([min(sortfewpoints(:,1)) max(sortfewpoints(:,1))]);

a.Color = 'b';
a.LineWidth = 1.5;
b.Color = 'r';
b.LineWidth = 1.5;
c.Color = 'b';
c.LineWidth = 1.5;
d.Color = 'r';
d.LineWidth = 1.5;
set(gca,'fontsize',14);
xlabel('% missing values that are MNAR');
ylabel('NS-kNN NRMSE - kNN NRMSE');

% Percentage at which kNN and NS-kNN error plots intersect. Uses InterX
% function by NS 2010 found on MathWorks File Exchange.

intersection = InterX([sorteverything(:,1)';sorteverything(:,2)'],[(0:1:100);zeros(1,101)]);
if isempty(intersection)
    0 % If no intesrection
else
    intersection(1,end)
end
