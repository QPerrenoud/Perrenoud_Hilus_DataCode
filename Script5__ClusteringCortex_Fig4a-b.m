clear, close all
chInputFolder   = 'D:\Perrenoud_Leclerc_etal\Analyses_Gyrus_Cortex\Perrenoud_Hilus_DataCode'; % < ---- CHANGE TO LOCAL PATH
cd(chInputFolder);
addpath(genpath('Utilities'))
chInputFile_1   = 'Input_Matrix_Hilus';
chInputFile_2   = 'Input_Matrix_Cortex300';
chOutputFolder  = 'Figures';
chSubFolder     = '5__ClusteringCortex_Fig4a-b';
% db1SizFig       = [100, 100, 1600, 900];
db1SizFig       = [50, 50, 1250, 600];

% Loads all the gyrus data
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
% % Removes VIP neurons
bl1VIP      = db2DataH0(:, 26) == 1;
db2DataH0(bl1VIP, :) = []; cNEURONS_H0(bl1VIP) = [];
% Removes non usefull parameters
in1IdxRem   = [7 15 22 26];
% in1IdxRem   = [7 15 22];
db2DataH0(:, in1IdxRem) = []; cLABELS_H0(in1IdxRem) = [];

% Seed random number generator
rng(1984)
%% Clusters Cortical GABAergic cells based on electrophysiological and molecular marker
db2ZDataH0      = zscore(db2DataH0);
db2Tree_Ctx     = linkage(db2ZDataH0, 'ward');
db1BIC_Ctx      = ClusterBIC(db2ZDataH0, db2Tree_Ctx, 10);
% db1BIC = ClusterLogLikelihood(db2ZDataH0, 10);
% [~, inNClu_Ctx] = min(db1BIC);
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
%% Clusters girus GABAergic cells based on electrophysiological and molecular marker
db2ZData        = zscore(db2Data);
db2Tree_Gyr     = linkage(db2ZData, 'ward');
inNClu_Ctx_Gyr      = 3;

% Calculates the centroid of ward clusters
in1CluWard_Gyr  = cluster(db2Tree_Gyr, 'maxclust', inNClu_Ctx_Gyr);
in1CluWard_Gyr  = SortClusterBySize(in1CluWard_Gyr);
db2CtrWard=zeros(inNClu_Ctx_Gyr, size(db2ZData, 2));
for iClu = 1:inNClu_Ctx_Gyr
    db2CtrWard(iClu,:) = mean(db2ZData(in1CluWard_Gyr == iClu, :));
end

% Does the K-Mean correction and prints the silhouette metrics
in1CluKM_Gyr    = kmeans(db2ZData, inNClu_Ctx_Gyr, 'start', db2CtrWard);

% If neeeded plots a dendrogram of the gyrus data for visual control
% figure,
% dendrogram(db2Tree_Gyr, size(db2ZData, 1), ...
%     'ColorThreshold', median(db2Tree_Gyr(end - inNClu_Ctx_Gyr + 1:end - inNClu_Ctx_Gyr + 2, 3)));

% Normalize the gyrus data to the cortex mean and s.d.
db1Mu_Ctx       = mean(db2DataH0);
db1SD_Ctx       = std(db2DataH0);
db2Data_ZCtx    = (db2Data - db1Mu_Ctx)./db1SD_Ctx;
db2Data_ZCtx(isnan(db2Data_ZCtx)) = 0;
%% Plots cortical aggregation dendrogram - and cortical and gyrus clusters on cortial principal components
close all
hFIG = []; iFig = 0; cFIG_LABEL = {};

% Initializes a color vector to plot clusters in different colors
db2Color_Ctx    = [.4 .8 0; 0 .5 0; .9 .45 0; .9 0 0];
cCOLOR_GYR      = {[255 222 23]./255, [57 181 74]./255, [185 82 159]./255};
db2Color_Gyr    = [];
for iCol = 1:length(cCOLOR_GYR); db2Color_Gyr = cat(1, db2Color_Gyr, cCOLOR_GYR{iCol}); end

% Plots the dendrogram for cortical neurons
iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
cFIG_LABEL = cat(1, cFIG_LABEL, ['Dendrogram_Cortex']);
subplot(4, 1, [1 2])
dendrogram(db2Tree_Ctx, size(db2ZDataH0, 1), 'ColorThreshold', median(db2Tree_Ctx(end - inNClu_Ctx + 1:end - inNClu_Ctx + 2, 3)));
subplot(4, 1, 3)
plot(2:10, db1BIC_Ctx(2:10), 'ko--'); ylabel('BIC'), xlabel('Number of clusters'), xlim([1 11]);
ylim([min(db1BIC_Ctx(2:10)) - .2 * range(db1BIC_Ctx(2:10)) max(db1BIC_Ctx(2:10)) + .2 * range(db1BIC_Ctx(2:10))]);
subplot(4, 1, 4)
plot(2:10, db1MinD_Ctx(2:10), 'ko--'); ylabel('Shortest Branch Length'), xlabel('Number of clusters'), xlim([1 11]);
ylim([min(db1MinD_Ctx(2:10)) - .2 * range(db1MinD_Ctx(2:10)) max(db1MinD_Ctx(2:10)) + .2 * range(db1MinD_Ctx(2:10))]);

% Plots the clusters on the two principal components
iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
cFIG_LABEL = cat(1, cFIG_LABEL, ['PCA_Cortex_Gyrus_AllPar']);
[db1Coeff, db1PC] = pca(db2ZDataH0);
subplot(1, 2, 1);
gscatter(db1PC(:,1), db1PC(:,2), in1CluKM_Ctx, db2Color_Ctx);
title('Cortical neurons in cortex space');
subplot(1, 2, 2);
gscatter(db2Data_ZCtx * db1Coeff(:, 1), db2Data_ZCtx * db1Coeff(:, 2), in1CluKM_Gyr, db2Color_Gyr);
title('Gyrus neurons in cortex space');
%% Plots Ephys and molecular parameters - appends values and statistical test to a text file
% Initialize figure and table variables

% Opens output file
hFID = fopen(fullfile(chOutputFolder, chSubFolder, 'parameters_cortex.txt'), 'w');

% Initializes colors for cortical neurons
cCOLOR_CTX  = {[.4 .8 0] [0 .5 0] [.9 .45 0] [.9 0 0]};

% Prints silhouette values
[inNObs, inPar] = size(db2ZDataH0);
db1Sil_1    = silhouette(db2ZDataH0, in1CluKM_Ctx);
fprintf(hFID, 'CORTICAL NEURONS : \r');
fprintf(hFID, '*************************\r\r');
fprintf(hFID, '-------- Silhouette -----------------------------\r');
fprintf(hFID, '%d cell out %d have a negative silhouette value\r', sum(db1Sil_1 < 0), inNObs);
fprintf(hFID, 'Mean silhouette: %.3f\r', mean(db1Sil_1));
fprintf(hFID, '\r');
fprintf(hFID, '\r');

% Plots values and significance for Ephys parameters
iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
cFIG_LABEL = cat(1, cFIG_LABEL, 'ParEphy_Cortex');
fprintf(hFID, '-------- Ephys ---------------------------------\r');
in1Phy      = 1:16; % ephys parameters indices
for iPar = in1Phy
    fprintf(hFID, '%s: ++++++++++++++++++++++++++++++++\r', cPAR{iPar});
    subplot(4, 4, iPar), hold on
    for iClu = 1:inNClu_Ctx
        bl1CluCells = in1CluKM_Ctx == iClu;
        inNCells    = sum(bl1CluCells);
        dbMean  = mean(db2DataH0(bl1CluCells, iPar));
        dbSEM   = sqrt(var(db2DataH0(bl1CluCells, iPar))/inNCells);
        errorbar(iClu, dbMean, dbSEM, dbSEM, 'o', 'Color', cCOLOR_CTX{iClu}, 'LineWidth', 2);
        fprintf(hFID, 'Cluster %d: %.3f (+- %.3f)\r', iClu, dbMean, dbSEM); % Print results to file
    end
    set(gca, 'XTick', 1:inNClu_Ctx);
    xlim([0 inNClu_Ctx + 1]);
    xlabel('Clusters');
    title(cPAR(iPar));
        
    % Prints the significance to test
    fprintf(hFID, '--------------Sigficance (Mann Witney):\r');
    for iClu = 1:inNClu_Ctx, fprintf(hFID, '\tC%d', iClu); end
    for iCl1 = 1:inNClu_Ctx - 1
        fprintf(hFID, '\r');
        fprintf(hFID, 'C%d', iCl1); for iClu = 1:iCl1, fprintf(hFID, '\t'); end
        for iCl2 = iCl1 + 1:inNClu_Ctx
            dbPVal = ranksum(db2DataH0(in1CluKM_Ctx == iCl1, iPar), db2DataH0(in1CluKM_Ctx == iCl2, iPar));
            fprintf(hFID, '\t%.3f', dbPVal);
        end
    end
    fprintf(hFID, '\r');
    fprintf(hFID, '\r');
end

in1Mol      = 17:24; % molecular parameters indices
% Plots the marker histograms of the clusters
iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
cFIG_LABEL  = cat(1, cFIG_LABEL, 'Histogram_Cortex');
cMARKER     = {'VGluT1', 'GAD', 'NOS1', 'PV', 'CR', 'NPY', 'SOM', 'CCK'};
bl1MrkIdx   = ismember(cLABELS_H0, cMARKER);
db2DatMark  = db2DataH0(:, bl1MrkIdx);


fprintf(hFID, '-------- Molecular markers -----------------------------\r');
fprintf(hFID, '\t');
for iMrk = 1:length(cMARKER), fprintf(hFID, '%s\t', cMARKER{iMrk}); end
fprintf(hFID, '\r');
[inRow, inCol] = FindPlotNum(inNClu_Ctx);
[db2Prop, in2NObs] = deal(nan(inNClu_Ctx, length(cMARKER)));
for iClu = 1:inNClu_Ctx
    subplot(inRow, inCol, iClu)
    db2Prop(iClu, :)    = mean(db2DatMark(in1CluKM_Ctx == iClu, :), 1);
    in2NObs(iClu, :)    = sum(in1CluKM_Ctx == iClu);
    bar(1:length(cMARKER), db2Prop(iClu, :))
    set(gca, 'XTick', 1:length(cMARKER), 'XTickLabel', cMARKER')
    ylim([0 1.2])
    title(sprintf('Cluster %d: N = %d', iClu, sum(in1CluKM_Ctx == iClu)));
    fprintf(hFID, 'Clu%d', iClu);
    for iMrk = 1:length(cMARKER), fprintf(hFID, '\t%.1f', db2Prop(iClu, iMrk) * 100); end
    fprintf(hFID, '\r');
end
fprintf(hFID, '\r');

% Prints the significance of molecular marker association
fprintf(hFID, '--------------Sigficance (Binomial Test):\r');
fprintf(hFID, '\t');
for iMrk = 1:length(cMARKER), fprintf(hFID, '%s\t', cMARKER{iMrk}); end
fprintf(hFID, '\r');
for iCl1 = 1:inNClu_Ctx - 1
    for iCl2 = iCl1 + 1:inNClu_Ctx
        db1PVal = BinomialTest2(in2NObs(iCl1, :), db2Prop(iCl1, :), in2NObs(iCl2, :), db2Prop(iCl2, :));
        fprintf(hFID, 'C%d-C%d', iCl1, iCl2);
        for iMrk = 1:length(cMARKER), fprintf(hFID, '\t%.3f', db1PVal(iMrk)); end
        fprintf(hFID, '\r');
    end
end
fprintf(hFID, '\r');
fprintf(hFID, '\r');

% Close the output writing file
fclose(hFID);
%% Makes sure that the exemple cells are in the proper cluster
cEXEMPLE_NEURON = {'SQ32', 'SQ101', 'SQ102'};
for iNeu = 1:length(cEXEMPLE_NEURON)
    in1NeuIdx = find(ismember(cNEURONS_H0, cEXEMPLE_NEURON(iNeu)));
    fprintf('%s\t: Cluster %d\r', cEXEMPLE_NEURON{iNeu}, in1CluKM_Ctx(in1NeuIdx));
end
%%
warning off
mkdir(chOutputFolder);
mkdir(fullfile(chOutputFolder, chSubFolder));
warning on
save(fullfile(chOutputFolder, chSubFolder, 'Workspace'), '-v7.3')
savefig(hFIG,fullfile(chOutputFolder, chSubFolder, 'Figures'));

for iFig = 1:length(hFIG)
    exportfig(hFIG(iFig) ,fullfile(chOutputFolder, chSubFolder, cFIG_LABEL{iFig}), 'Color', 'rgb', 'Renderer', 'painters');
    set(hFIG(iFig), 'Units', 'centimeters')
    db1Pos = get(hFIG(iFig), 'Position');
    set(hFIG(iFig), 'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', 'PaperSize', [db1Pos(3) db1Pos(4)])
    print(hFIG(iFig), fullfile(chOutputFolder, chSubFolder, cFIG_LABEL{iFig}), '-dpdf')
end