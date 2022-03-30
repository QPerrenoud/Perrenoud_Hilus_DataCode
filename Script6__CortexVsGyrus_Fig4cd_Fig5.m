clear, close all
chInputFolder   = 'D:\Perrenoud_Leclerc_etal\Perrenoud_Hilus_DataCode'; % < ---- CHANGE TO LOCAL PATH
cd(chInputFolder);
addpath(genpath('Utilities'))
chInputFile_1   = 'Input_Matrix_Hilus';
chInputFile_2   = 'Input_Matrix_Cortex300';
chOutputFolder  = 'Figures';
chSubFolder     = '6__CortexVsGyrus_Fig4cd_Fig5';
% db1SizFig       = [100, 100, 1600, 900];
db1SizFig       = [50, 50, 1250, 600];

% Loads all the data
[~, cNEURONS]   = xlsread(fullfile(chInputFolder, chInputFile_1), 'A2:A149');
[~, cLABELS]    = xlsread(fullfile(chInputFolder, chInputFile_1), 'K1:CN1');
db2DataAll      = xlsread(fullfile(chInputFolder, chInputFile_1), 'K2:CN149');
% for iPar = 1:length(cLABELS); fprintf('%d : %s\r', iPar, cLABELS{iPar}); end

% Construct the parameter matrix
in1Idx      = [1:4 7 9 13:16 19:20 24 28 39 40 41:42 44 46:48 50:51]; %Loads the parameters
cPAR        = cLABELS(in1Idx); %
db2Data     = db2DataAll(:, in1Idx);
% Merges the expression of the two gads
cPAR{18}    = 'GAD';
db2Data(:, 18) = db2DataAll(:,42) | db2DataAll(:, 43);
% Selects GAD neurons
bl1GAD      = db2Data(:, 18) == 1;
db2Data     = db2Data(bl1GAD, :); cNEURONS = cNEURONS(bl1GAD);
% Loads the cortex data
[~, cNEURONS_H0]    = xlsread(fullfile(chInputFolder, chInputFile_2), 2, 'A2:A301');
[~, cLABELS_H0]     = xlsread(fullfile(chInputFolder, chInputFile_2), 2, 'B1:AE1');
db2DataH0           = xlsread(fullfile(chInputFolder, chInputFile_2), 2, 'B2:AE301');
% Removes VIP neurons
bl1VIP      = db2DataH0(:, 26) == 1;
db2DataH0(bl1VIP, :) = []; cNEURONS_H0(bl1VIP) = [];
% Removes non usefull parameters
in1IdxRem   = [7 15 22 26];
db2DataH0(:, in1IdxRem) = [];

% Seed random number generator
rng(1984);
%% Clusters Cortical GABAergic cells based on electrophysiological and molecular marker
db2ZDataH0      = zscore(db2DataH0);
db2Tree_Ctx     = linkage(db2ZDataH0, 'ward');
db1BIC_Ctx      = ClusterBIC(db2ZDataH0, db2Tree_Ctx, 10);
% db1BIC = ClusterLogLikelihood(db2ZDataH0, 10);
% [~, inNClu] = min(db1BIC);
db1MinD_Ctx     = ClusterDist(db2Tree_Ctx, 10);
inNClu_Ctx      = 4;

% Calculates the centroid of ward clusters
in1CluWard_Ctx  = cluster(db2Tree_Ctx, 'maxclust', inNClu_Ctx);
db2CtrWard=zeros(inNClu_Ctx, size(db2ZDataH0, 2));
for iClu = 1:inNClu_Ctx
    db2CtrWard(iClu,:) = mean(db2ZDataH0(in1CluWard_Ctx == iClu, :));
end

% Does the K-Mean correction and prints the silhouette metrics
in1CluKM_Ctx    = kmeans(db2ZDataH0, inNClu_Ctx, 'start', db2CtrWard);
%% Clusters Cortical GABAergic cells based on electrophysiological and molecular marker
db2ZData        = zscore(db2Data);
db2Tree_Gyr     = linkage(db2ZData, 'ward');
inNClu_Gyr      = 3;

% Calculates the centroid of ward clusters
in1CluWard_Gyr  = cluster(db2Tree_Gyr, 'maxclust', inNClu_Gyr);
in1CluWard_Gyr  = SortClusterBySize(in1CluWard_Gyr);
db2CtrWard=zeros(inNClu_Gyr, size(db2ZData, 2));
for iClu = 1:inNClu_Gyr
    db2CtrWard(iClu,:) = mean(db2ZData(in1CluWard_Gyr == iClu, :));
end

% Does the K-Mean correction and prints the silhouette metrics
in1CluKM_Gyr    = kmeans(db2ZData, inNClu_Gyr, 'start', db2CtrWard);

% If neeeded plots a dendrogram of the gyrus data for visual control
% figure,
% dendrogram(db2Tree_Gyr, size(db2ZData, 1), ...
%     'ColorThreshold', median(db2Tree_Gyr(end - inNClu_Gyr + 1:end - inNClu_Gyr + 2, 3)));
%% Test the angle between the first PC of different substsets of Hilar and Cortical GABA neurons
cPAR_GP     = {'All', 'NoRMP_NoCSat'};
cREM_PAR    = {[], [1 9],};
inNPGp      = length(cPAR_GP);
cIN_GP      = {'GAD' 'PV' 'SOM' 'NPY' 'NOS1', 'CCK'};
cCRIT_IDX   = {18, [18 20], [18 23] [18 22], [18 19], [18 24]};
inNIGp      = length(cIN_GP);

% Initializes a color vector to plot clusters in different colors
    % For Marker expression on scatter plots
in1MrkIdx_All   = [20 23 22 19];
db2Color_Mkr    = [1 0 0; 0 1 0; 1 .5 0; .5 .5 1];
    % For Cluster allocation on scatter plots
db2Color_Ctx    = [.4 .8 0; 0 .5 0; .9 .45 0; .9 0 0];
cCOLOR_GYR      = {[255 222 23]./255, [57 181 74]./255, [185 82 159]./255};
db2Color_Gyr    = [];
for iCol = 1:length(cCOLOR_GYR); db2Color_Gyr = cat(1, db2Color_Gyr, cCOLOR_GYR{iCol}); end
    % For Quantiles
cQ_Color = {[.7 .7 .7], [.4 .4 .4], [.4 .4 .4], [.7 .7 .7]};

for iPGp = 1:inNPGp
    %Makes the figure output folder
    warning off
    %mkdir(chOutputFolder);
    %mkdir(fullfile(chOutputFolder, chSubFolder));
    %mkdir(fullfile(chOutputFolder, chSubFolder, cPAR_GP{iPGp}));
    chGroupFolder = fullfile(chOutputFolder, chSubFolder, cPAR_GP{iPGp});
    mkdir(fullfile(chGroupFolder));
    
    %Opens a file to print the stats
    hFID = fopen(fullfile(chGroupFolder, 'Stats.txt'), 'w');
    
    clear hFIG cFIG_LABEL db2H db1H
    iFig = 0;
    
    %Initializes aggregation variables
    [db1DotP, db1PVal_DotP, db1CMahalDist, db1PVal_Dist] = deal(nan(1, inNIGp));
    [db2CI_DotP, db2CI_Dist] = deal(nan(5, inNIGp));
    
    %Communicates
    fprintf('%s parameters:\n---------------------------------------------\r', ...
        cPAR_GP{iPGp})
    fprintf(hFID, '%s parameters:\n---------------------------------------------\r', ...
        cPAR_GP{iPGp});
    
    %Loops over marker of interest
    for iIGp = 1:inNIGp
        % Select the data and performs the test
        bl1SelH0        = all(db2DataH0(:, cCRIT_IDX{iIGp}), 2);
        db2DataH0_Tst   = db2DataH0(bl1SelH0, :);
        db2DataH0_Tst(:, [cCRIT_IDX{iIGp} cREM_PAR{iPGp}]) = [];
        bl1SelTst       = all(db2Data(:, cCRIT_IDX{iIGp}), 2);
        db2Data_Tst     = db2Data(bl1SelTst, :);
        db2Data_Tst(:, [cCRIT_IDX{iIGp} cREM_PAR{iPGp}]) = [];
        [db1DotP(iIGp), db1PVal_DotP(iIGp), db2CI_DotP(:, iIGp), db2Score_H0, ...
            db2Score_Tst, db1VarExp_H0, db1VarExp_Test] = PCAAngleTest(db2DataH0_Tst, db2Data_Tst);
        [db1CMahalDist(iIGp), db1PVal_Dist(iIGp), db2CI_Dist(:, iIGp)] = ...
            CentroidMahalDistTest(db2DataH0_Tst, db2Data_Tst);
        
        %Communicates
        fprintf('%s\tPC1 DotProduct\t\t\t\t%.3f (p = %.4f)\n\tCentroid Mahalanobis Dist\t%.2f (p = %.4f)\n', ...
            cIN_GP{iIGp}, db1DotP(iPGp), db1PVal_DotP(iIGp), db1CMahalDist(iIGp), db1PVal_Dist(iIGp));
        fprintf('Variance Explained H0:\t%.2f\nVariance Explained Test:\t%.2f\n', ...
            db1VarExp_H0(1) * 100, db1VarExp_Test(1) * 100);
        fprintf(hFID, '%s\tPC1 DotProduct\t\t\t\t%.3f (p = %.4f)\n\tCentroid Mahalanobis Dist\t%.2f (p = %.4f)\n', ...
            cIN_GP{iIGp}, db1DotP(iPGp), db1PVal_DotP(iIGp), db1CMahalDist(iIGp), db1PVal_Dist(iIGp));
        fprintf(hFID, 'Variance Explained H0:\t%.2f\nVariance Explained Test:\t%.2f\n', ...
            db1VarExp_H0(1) * 100, db1VarExp_Test(1) * 100);
%         continue
        
        % Plots the figure
        cMRK_COL = {[1 0 0], [0 1 0], [1 .5 0]};
        iFig = iFig + 1;
        hFIG(iFig) = figure('Position', db1SizFig);
        cFIG_LABEL{iFig} = cIN_GP{iIGp};
        in1MrkIdx   = in1MrkIdx_All;
        db2Color    = db2Color_Mkr;
        %bl1Keep = ~ismember(in1MrkIdx_All, cCRIT_IDX{iIGp});
        %in1MrkIdx   = in1MrkIdx_All(bl1Keep);
        %db2Color    = db2Color_Mkr(bl1Keep, :);
        
        subplot(2, 2, 1); pause(.1);
        PieScatter(db2Score_H0(:, 1), db2Score_H0(:, 2), db2DataH0(bl1SelH0, in1MrkIdx), db2Color)
        xlabel(sprintf('PC1 (%.2fpct)', db1VarExp_H0(1) * 100)), ylabel(sprintf('PC2 (%.2fpct)', db1VarExp_H0(2) * 100)), 
        title(sprintf('%s : Cortex IN in Cortex Space (n = %d)', cIN_GP{iIGp}, sum(bl1SelH0)));
        db1XL = xlim; db1YL = ylim;
        subplot(2, 2, 2); pause(.1);
        PieScatter(db2Score_Tst(:, 1), db2Score_Tst(:, 2), db2Data(bl1SelTst, in1MrkIdx), db2Color, db1XL, db1YL)
        xlabel(sprintf('PC1 (%.2fpct)', db1VarExp_H0(1) * 100)), ylabel(sprintf('PC2 (%.2fpct)', db1VarExp_H0(2) * 100)),
        title(sprintf('Gyrus INs in Cortex Space (n = %d, DotProduct : %.2f, p = %.3f; Mhl Dist : %.3f , p = %.3f;)',...
            sum(bl1SelTst), db1DotP(iIGp), db1PVal_DotP(iIGp), db1CMahalDist(iIGp), db1PVal_Dist(iIGp)));
        subplot(2, 2, 3)
        xlabel('PC1'), ylabel('PC2'),
            % calculates the CC line
        dbX_C_H0    = mean(db2Score_H0(:, 1));  dbY_C_H0    = mean(db2Score_H0(:, 2));
        dbX_C_Tst   = mean(db2Score_Tst(:, 1)); dbY_C_Tst   = mean(db2Score_Tst(:, 2));
        dbA         = (dbY_C_H0 - dbY_C_Tst) ./ (dbX_C_H0 - dbX_C_Tst);
        dbB         = dbY_C_H0 - (dbX_C_H0 * dbA);
        fLINE       = @(x) dbA * x + dbB;
        gscatter(db2Score_H0(:, 1), db2Score_H0(:, 2), in1CluKM_Ctx(bl1SelH0), db2Color_Ctx); hold on 
        plot(dbX_C_H0, dbY_C_H0, 'kx'); % adds the centroid
        plot(dbX_C_Tst, dbY_C_Tst, 'kx'); 
        plot(db1XL, fLINE(db1XL), 'k--'); % adds the axe connecting the centroids
        xlim(db1XL), ylim(db1YL);
        subplot(2, 2, 4);
        gscatter(db2Score_Tst(:, 1), db2Score_Tst(:, 2), in1CluKM_Gyr(bl1SelTst), db2Color_Gyr); hold on
        plot(dbX_C_Tst, dbY_C_Tst, 'kx'); % adds the centroid
        plot(dbX_C_H0, dbY_C_H0, 'kx'); 
        plot(db1XL, fLINE(db1XL), 'k--'); % adds the axe connecting the centroids
        xlim(db1XL), ylim(db1YL);
        
        %Plots a histogram of the points projected on the axis connecting
        %their centroids
        [db1Proj_H0, db1Proj_Test] = CCAxeProjection(db2DataH0_Tst, db2Data_Tst);
        db1Bin = -1.5:.25:2.5;
        db2H_Ctx = []; cLEGEND_CTX = {};
        for iClu = 1:inNClu_Ctx 
            db1H = hist(db1Proj_H0(in1CluKM_Ctx(bl1SelH0) == iClu), db1Bin);
            db2H_Ctx = cat(1, db2H_Ctx, db1H);
            cLEGEND_CTX = cat(2, cLEGEND_CTX, {sprintf('Ctx%d', iClu)});
        end
        db2H_Gyr = []; cLEGEND_GYR = {};
        for iClu = 1:inNClu_Gyr
            db1H = hist(db1Proj_Test(in1CluKM_Gyr(bl1SelTst) == iClu), db1Bin);
            db2H_Gyr = cat(1, db2H_Gyr, db1H);
            cLEGEND_GYR = cat(2, cLEGEND_GYR, {sprintf('Gyr%d', iClu)});
        end
        iFig = iFig + 1;
        hFIG(iFig) = figure('Position', db1SizFig);
        cFIG_LABEL{iFig} = [cIN_GP{iIGp} '_CCProj'];
            % Ctx and Gyrus
        hAX(1) = subplot(1, 3, 1);
        hB = bar(db1Bin, [db2H_Ctx; db2H_Gyr]', 'stacked');
        db2H_Color = [db2Color_Ctx; db2Color_Gyr];
        for iClu = 1:length(hB); set(hB(iClu), 'FaceColor', db2H_Color(iClu, :)); end
        legend(hB, cat(2, cLEGEND_CTX, cLEGEND_GYR));
        set(gca, 'XTick', [0 1], 'XTickLabel', {'C-Ctx', 'C-Gyr'})
        title('C-C Projection')
            % Ctx
        hAX(2) = subplot(1, 3, 2);
        hB = bar(db1Bin, db2H_Ctx', 'stacked');
        db2H_Color = [db2Color_Ctx];
        for iClu = 1:length(hB); set(hB(iClu), 'FaceColor', db2H_Color(iClu, :)); end
        legend(hB, cLEGEND_CTX);
        set(gca, 'XTick', [0 1], 'XTickLabel', {'C-Ctx', 'C-Gyr'})
        title('C-C Projection')
            % Gyr
        hAX(3) = subplot(1, 3, 3);
        hB = bar(db1Bin, db2H_Gyr', 'stacked');
        db2H_Color = [db2Color_Gyr];
        for iClu = 1:length(hB); set(hB(iClu), 'FaceColor', db2H_Color(iClu, :)); end
        legend(hB, cLEGEND_GYR);
        set(gca, 'XTick', [0 1], 'XTickLabel', {'C-Ctx', 'C-Gyr'})
        title('C-C Projection')
        linkaxes(hAX, 'y');
        
        %Test the statistical difference of all parameters between the two data
        %sets and order them
        cPAR_FIG = cPAR; cPAR_FIG([cCRIT_IDX{iIGp} cREM_PAR{iPGp}]) = [];
        db1P_RkSm = nan(1, length(cPAR_FIG));
        for iPar = 1:length(cPAR_FIG)
            db1P_RkSm(iPar) = ranksum(db2DataH0_Tst(:, iPar), db2Data_Tst(:, iPar));
        end
        [db1SortP, in1Idx] = sort(db1P_RkSm);
        
        iFig = iFig + 1;
        hFIG(iFig) = figure('Position', db1SizFig);
        cFIG_LABEL{iFig} = [cIN_GP{iIGp} '_SigMarker'];
        plot(1:length(cPAR_FIG), db1SortP, 'ko')
        set(gca, 'XTick', 1:length(cPAR_FIG), 'XTickLabel', cPAR_FIG(in1Idx), 'YScale', 'log')
        ylim([min(db1SortP)./10 1])
        ylabel('PValue (Ranksum Test)')
        title(sprintf('%d significant difference', sum(db1SortP < .05)));
        
        %Generates 2d plots of the most significant difference
        iFig = iFig + 1;
        hFIG(iFig) = figure('Position', db1SizFig);
        cFIG_LABEL{iFig} = [cIN_GP{iIGp} '_SigMarkerComp'];
        for iPlt = 1:12
            subplot(3, 4, iPlt)
            inIdxPlt    = in1Idx(iPlt); hold on
            db1DatH0    = db2DataH0_Tst(:, inIdxPlt);
            db1Dat      = db2Data_Tst(:, inIdxPlt);
            db1Mean     = [mean(db1DatH0) mean(db1Dat)];
            db1Err      = [std(db1DatH0) std(db1Dat)];
            plot(1, db1DatH0, 'x', 'Color', [.5 .5 .5]);
            plot(2, db1Dat, 'x', 'Color', [.5 .5 .5]);
            errorbar([1 2], db1Mean, db1Err, db1Err, 'ko', 'LineWidth', 2)
            set(gca, 'XTick', [1 2], 'XTickLabel', {'Cortex', 'Gyrus'});
            xlim([0 3]);
            %ylabel(cPAR_FIG{inIdxPlt})
            title(sprintf('%s : p = %.4f', cPAR_FIG{inIdxPlt}, db1SortP(iPlt)));
            fprintf('\n\n%s : p = %.4f', cPAR_FIG{inIdxPlt}, db1SortP(iPlt));
            fprintf('\nCortex (mean +/- S.D):\t%.2f +/- %.2f', db1Mean(1), db1Err(1));
            fprintf('\nHilus (mean +/- S.D):\t%.2f +/- %.2f', db1Mean(2), db1Err(2));
            fprintf(hFID, '\n\n%s : p = %.4f', cPAR_FIG{inIdxPlt}, db1SortP(iPlt));
            fprintf(hFID, '\nCortex (mean +/- S.D):\t%.2f +/- %.2f', db1Mean(1), db1Err(1));
            fprintf(hFID, '\nHilus (mean +/- S.D):\t%.2f +/- %.2f', db1Mean(2), db1Err(2));
        end
        
        pause(.1)
        SaveFig(chGroupFolder, hFIG, cFIG_LABEL);
        close all
        clear hFIG cFIG_LABEL
        iFig = 0;
    end
    fclose(hFID);
    fprintf('\n');
    
    %Plots a figure of the dot products of the first PCs and their p-value
    iFig = iFig + 1;
    hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL{iFig} = 'Dot Product All';
        %Mahalanobis distance
    subplot(1, 2, 1)
    plot(1:inNIGp, db1DotP, 'ko'); hold on
    for iIGp = 1:inNIGp
        db1Y = [-.25 .25 .25 -.25] + iIGp;
        for iQ = 1:size(db2CI_DotP) - 1
            fill(db1Y, [[1 1] * db2CI_DotP(iQ, iIGp) [1 1] * db2CI_DotP(iQ + 1, iIGp)] , ...
                cQ_Color{iQ}, 'LineStyle', 'none', 'FaceAlpha', .3)
        end
        plot([-.25 .25] + iIGp, [1 1] * db2CI_DotP(3, iIGp), 'r');
    end
    xlim([0 inNIGp + 1])
    set(gca, 'XTick', 1:inNIGp, 'XTickLabel', cIN_GP);
    xlabel('Marker Group')
    title('Dot Product of PCA angle');
        %P-Value
    subplot(1, 2, 2)
    plot(1:inNIGp, db1PVal_DotP, 'ko');
    xlim([0 inNIGp + 1])
    ylim([0 max(db1PVal_DotP) + .1])
    set(gca, 'XTick', 1:inNIGp, 'XTickLabel', cIN_GP);
    xlabel('Marker Group')
    title('P-Value of Dot Product');
    
    %Plots a figure of the mahalanobis distances and the p-value of
    %belonging to the H0 group for each interneuron group
    iFig = iFig + 1;
    hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL{iFig} = 'Mahalanobis Distance All';
        %Mahalanobis distance
    subplot(1, 2, 1)
    plot(1:inNIGp, db1CMahalDist, 'ko'); hold on
    for iIGp = 1:inNIGp
        db1Y = [-.25 .25 .25 -.25] + iIGp;
        for iQ = 1:size(db2CI_Dist) - 1
            fill(db1Y, [[1 1] * db2CI_Dist(iQ, iIGp) [1 1] * db2CI_Dist(iQ + 1, iIGp)] , ...
                cQ_Color{iQ}, 'LineStyle', 'none', 'FaceAlpha', .3)
        end
        plot([-.25 .25] + iIGp, [1 1] * db2CI_Dist(3, iIGp), 'r');
    end
    xlim([0 inNIGp + 1])
    set(gca, 'XTick', 1:inNIGp, 'XTickLabel', cIN_GP);
    xlabel('Marker Group')
    title('Mahalanobis Distance to Cortex Centroid');
        %P-Value
    subplot(1, 2, 2)
    plot(1:inNIGp, db1PVal_Dist, 'ko');
    xlim([0 inNIGp + 1])
    ylim([0 max(db1PVal_Dist) + .1])
    set(gca, 'XTick', 1:inNIGp, 'XTickLabel', cIN_GP);
    xlabel('Marker Group')
    title('P-Value of Centroid Distance');
    
    SaveFig(chGroupFolder, hFIG, cFIG_LABEL);
    close all
    warning on
end
