clear, close all
chInputFolder   = 'C:\Users\perre\Dropbox\Leclerc et al. 2018\Brain Structure Function\New Analysis';
% chInputFolder   = '/home/quentin/Dropbox/Leclerc et al. 2018/Brain Structure Function/New Analysis';
cd(chInputFolder);
addpath(genpath('Utilities'))
chInputFile     = 'Input_Matrix_NewClustering';
chOutputFolder  = 'Figures';
chSubFolder     = 'GAD_ParameterSets';
% db1SizFig       = [100, 100, 1600, 900];
db1SizFig       = [50, 50, 1250, 600];

% Loads all the data
[~, cNEURONS]   = xlsread(fullfile(chInputFolder, chInputFile), 'A2:A149');
[~, cLABELS]    = xlsread(fullfile(chInputFolder, chInputFile), 'K1:CN1');
db2DataAll      = xlsread(fullfile(chInputFolder, chInputFile), 'K2:CN149');
% for iPar = 1:length(cLABELS); fprintf('%d : %s\r', iPar, cLABELS{iPar}); end

% Construct the parameter matrix
in1Idx      = [1:4 7 9 13:16 19:20 24 28 39 40 41:42 44 46:48 50:51]; %Loads the parameters
cPAR        = cLABELS(in1Idx); %
db2Data     = db2DataAll(:, in1Idx);
cPAR{18}    = 'GAD';
db2Data(:, 18) = db2DataAll(:,42) | db2DataAll(:, 43); % Merges the expression of the two gads
db2ZData    = zscore(db2Data);

% Seed random number generator
rng(1984)
%% Clusters GABAergic cells based on electrophysiological or molecular marker
% Defines sub populations of interest for the analysis
cCLU_POP    = {'GAD'};
in1CluNum   = 2;
cSEL        = {db2Data(:, 18) == 1};
iPop        = 1;

% Creates a boolean vector indicating if parameters are molecular
in1Phy      = 1:16;
in1Mol      = 17:24;
in1Sig      = [1:6 11 13 15:17 21:23];
in1SigPhy   = [1:6 11 13 15:16];
in1SigMol   = [17 21:23];
cPAR_SEL    = {true(size(cPAR)), ismember(cPAR, cPAR(in1Phy)), ismember(cPAR, cPAR(in1Mol)), ...
    ismember(cPAR, cPAR(in1Sig)), ismember(cPAR, cPAR(in1SigPhy)), ismember(cPAR, cPAR(in1SigMol))};
cCLU_MSEL   = {'All', 'Ephy', 'Mol', 'Sig', 'SigPhy', 'SigMol'};
inNSel      = length(cPAR_SEL);

% Initializes a color vector to plot clusters in different colors
cCOLOR = {[255 222 23]./255, [57 181 74]./255, [185 82 159]./255, [.4 .8 0], [0 .5 0], [.9 .45 0], [.9 0 0], ...
    [.9 0 .45], [0 0 .9], [0 .45 .9], [.45 0 .9], [.45 .9 0], [0 .9 .45]};
% rng(1984);
for iItr = 1:10
    cCOLOR = cat(2, cCOLOR, {rand(1, 3)});
end

% % Initialize figure and table variables
hFIG = []; iFig = 0; cFIG_LABEL = {}; inNPltCol = 3;
% cTAB = {}; iTab = 0; cTAB_LABEL = {};

% Intiaalize the output variable
in2_CluKMean = nan(sum(cSEL{iPop}), inNSel);

% Loops through condition
for iPSl = 1:inNSel
    % Selects the fraction of the data
%     db2ZData_PSel       = db2ZData(cSEL{iPop}, cPAR_SEL{iPSl}); %Out 2020-02-25 
    db2Data_PSel        = db2Data(cSEL{iPop}, cPAR_SEL{iPSl});
    db2ZData_PSel       = zscore(db2Data_PSel);
    db2DataAll_PSel     = db2DataAll(cSEL{iPop}, :);
    [inNObs, inNVar]    = size(db2ZData_PSel);
    
    % Clusters cell according to Ephy or molecular markers
    % Clusters the data with ward methods
    db2Tree     = linkage(db2ZData_PSel, 'ward');
%     db2Inconst  = inconsistent(db2Tree); % Use the inconsistency to determine the number of cluster
%     db1Inconst  = db2Inconst(end:-1:1, 4); %Inverts the order of the rows
%     inNClu      = max(find(db1Inconst < 1, 1, 'first'), 2); % Find the first node with inconsistency inferior to one
%     inNClu      = 2;
%     [~, inNClu] = min(ClusterBIC(db2ZData_PSel, db2Tree, 10));
    [~, inNClu] = max(ClusterDist(db2Tree, 10));

    % Plots the dendrogram
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Dendrogram_' cCLU_MSEL{iPSl}]);
    if inNClu == 1
        dendrogram(db2Tree, size(db2ZData_PSel, 1));
    else
        dendrogram(db2Tree, size(db2ZData_PSel, 1), ...
            'ColorThreshold', median(db2Tree(end - inNClu + 1:end - inNClu + 2, 3)));
    end
    title(cCLU_MSEL{iPSl});
    
    if inNClu == 1; continue; end
    
    % Calculates the centroid of ward clusters
    in1CluWard  = cluster(db2Tree, 'maxclust', inNClu);
    % Ensures that the clusters are ordered by number of cells
    in1CluWard = SortClusterBySize(in1CluWard);   
    % Initializes centroids
    db2CtrWard=zeros(inNClu, size(db2ZData_PSel, 2));
    for iClu = 1:inNClu
        db2CtrWard(iClu,:) = mean(db2ZData_PSel(in1CluWard == iClu, :));
    end
    
    % Does the K-Mean correction and prints the silhouette metrics
    in1CluKM    = kmeans(db2ZData_PSel, inNClu, 'start', db2CtrWard);
    db1Sil      = silhouette(db2ZData_PSel, in1CluKM);
    fprintf('%d cellules sur %d ont une valeur de silhouette negative\r', length(db1Sil(db1Sil < 0)), inNObs);
    fprintf('Silhouette moyenne: %.3f\r', mean(db1Sil));
    
    %Plots the silhouette
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Silhouette_' cCLU_MSEL{iPSl}]);
    silhouette(db2ZData_PSel, in1CluKM);
    
    % Plost the value of electrophysiological parameters
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['ParEphy_' cCLU_MSEL{iPSl}]);
    for iPar = in1Phy
        subplot(4, 4, iPar), hold on
        for iClu = 1:inNClu
            bl1CluCells = in1CluKM == iClu;
            inNCells    = sum(bl1CluCells);
            dbMean  = mean(db2Data(bl1CluCells, iPar));
            dbSEM   = sqrt(var(db2Data(bl1CluCells, iPar))/inNCells);
            errorbar(iClu, dbMean, dbSEM, dbSEM, 'o', 'Color', cCOLOR{iClu}, 'LineWidth', 2);
        end
        set(gca, 'XTick', 1:inNClu);
        xlim([0 inNClu + 1]);
        xlabel('Clusters');
        title(cPAR(iPar));
    end
    
    % Plots the marker histograms of the clusters
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Histogram_' cCLU_MSEL{iPSl}]);
    cMARKER     = {'VGluT1', 'GAD65', 'GAD67', 'NOS1', 'CA', 'PV', 'CR', 'NPY', 'VIP', 'SOM', 'CCK'};
    in1MrkIdx   = cellfun(@(x) find(ismember(cLABELS, x)), cMARKER);
    db2DatMark  = db2DataAll_PSel(:, in1MrkIdx);
    db2DatMark(:, 2) = db2DatMark(:, 2) | db2DatMark(:, 3);
    cMARKER(3) = []; cMARKER{2} = 'GAD';
    db2DatMark(:, 3) = [];
    if inNClu <= inNPltCol, inRow = 1; inCol = inNClu;
    else, inRow = ceil(inNClu/inNPltCol); inCol = inNPltCol; end
    for iClu = 1:inNClu
        subplot(inRow, inCol, iClu)
        bar(1:length(cMARKER), mean(db2DatMark(in1CluKM == iClu, :)))
        set(gca, 'XTick', 1:length(cMARKER), 'XTickLabel', cMARKER')
        ylim([0 1.2])
        title(sprintf('Cluster %d: N = %d', iClu, sum(in1CluKM == iClu)))
    end
    
    % Calculates the centroid of ward and the k-mean for two clusters to
    % compare the results of the different parameter sets
    inNClu = 3;
    in1CluWard  = cluster(db2Tree, 'maxclust', inNClu);
    in1CluWard  = SortClusterBySize(in1CluWard);
    db2CtrWard=zeros(inNClu, size(db2ZData_PSel, 2));
    for iClu = 1:inNClu
        db2CtrWard(iClu,:) = mean(db2ZData_PSel(in1CluWard == iClu, :));
    end
    in1_3CKM = kmeans(db2ZData_PSel, inNClu, 'start', db2CtrWard);
    
    % Rearrange clusters so that they correspond to one another across
    % parameter change
    if iPSl == 1, in1_KM_Srt = in1_3CKM; %If PSel == 1, takes the clustering as a reference
    else % finds which of the reference clusters has the most cells in the present cluster
        in1_KM_Srt = nan(size(in1_3CKM));
        for iClu = 1:inNClu
            in1Idx = in2_CluKMean(:, 1) == iClu;
            db1Hist = hist(in1_3CKM(in1Idx), 1:inNClu);
            [~, inMaxClu] = max(db1Hist);
            in1_KM_Srt(in1_3CKM == inMaxClu) = iClu;
        end
    end
    in2_CluKMean(:, iPSl) = in1_KM_Srt;
    
%     pause
    pause(0.1), %close all
end

% Sorts the cluster allocation and plots
[~, in1ObSrt]   = sort(in2_CluKMean(:, 1));
in2_2CKM_Srt    = in2_CluKMean(in1ObSrt, :); 
% Creates a colormap
db2ColorMap = [];
for iCol = 1:inNClu, db2ColorMap = cat(1, db2ColorMap, cCOLOR{iCol}); end
iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
cFIG_LABEL = cat(1, cFIG_LABEL, 'Cluster_Comparision');
imagesc(1:inNSel, 1:inNObs, in2_2CKM_Srt), colormap(db2ColorMap);
set(gca, 'XTick', 1:inNSel, 'XTickLabel', cCLU_MSEL);
ylabel('Neuron #');
title('GAD neurons : Cluster allocation as a function of parameter set')
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