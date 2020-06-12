clear, close all
chInputFolder   = 'C:\Users\perre\Dropbox\Leclerc et al. 2018\Brain Structure Function\New Analysis';
% chInputFolder   = '/home/quentin/Dropbox/Leclerc et al. 2018/Brain Structure Function/New Analysis';
% chInputFolder   = '/Users/Quentin/Dropbox/Leclerc et al. 2018/Brain Structure Function/New Analysis';
cd(chInputFolder);
addpath(genpath('Utilities'))
chInputFile     = 'Input_Matrix_NewClustering';
chOutputFolder  = 'Figures';
chSubFolder     = 'C2_EphySubclustering';
% db1SizFig       = [100, 100, 1600, 900];
db1SizFig       = [50, 50, 1250, 600];

% Loads all the data
[~, cNEURONS]   = xlsread(fullfile(chInputFolder, chInputFile), 'A2:A149');
[~, cLABELS]    = xlsread(fullfile(chInputFolder, chInputFile), 'K1:CN1');
db2DataAll      = xlsread(fullfile(chInputFolder, chInputFile), 'K2:CN149');
% for iPar = 1:length(cLABELS); fprintf('%d : %s\r', iPar, cLABELS{iPar}); end

% Construct the parameter matrixsort(sum(in2_2CluKMean, 2))
in1Idx      = [1:4 7 9 13:16 19:20 24 28 39 40 41:42 44 46:48 50:51]; %Loads the parameters
cPAR        = cLABELS(in1Idx); %
db2Data     = db2DataAll(:, in1Idx);
cPAR{18}    = 'GAD';
db2Data(:, 18) = db2DataAll(:,42) | db2DataAll(:, 43); % Merges the expression of the two gads

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
in1Mol      = 17:24;% Print the list of the cells in each cluster
in1Sig      = [1 2 5 6 15 16 22 23 24];
in1SigPhy   = [1 2 5 6 15 16];
in1SigMol   = [22 23 24];
cPAR_SEL    = {true(size(cPAR)), ismember(cPAR, cPAR(in1Phy)), ismember(cPAR, cPAR(in1Mol))};
cCLU_MSEL   = {'All', 'Ephy', 'Mol'};
inNSel      = length(cPAR_SEL);

% Initializes a color vector to plot clusters in different colors
cCOLOR = {[255 222 23]./255, [57 181 74]./255, [185 82 159]./255, [.4 .8 0], [0 .5 0], [.9 .45 0], [.9 0 0], ...
    [.9 0 .45], [0 0 .9], [0 .45 .9], [.45 0 .9], [.45 .9 0], [0 .9 .45]};
% rng(1984);
for iItr = 1:10
    cCOLOR = cat(2, cCOLOR, {rand(1, 3)});
end
%% Clusters the data
% Selects the fraction of the data
db2Data_GAD     = db2Data(cSEL{iPop}, :);
db2ZData_GAD    = zscore(db2Data_GAD);

% Clusters cell according to Ephy or molecular markers
% Clusters the data with ward methods
db2Tree     = linkage(db2ZData_GAD, 'ward');
%     db2Inconst  = inconsistent(db2Tree); % Use the inconsistency to determine the number of cluster
%     db1Inconst  = db2Inconst(end:-1:1, 4); %Inverts the order of the rows
%     inNClu      = max(find(db1Inconst < 1, 1, 'first'), 2); % Find the first node with inconsistency inferior to one
inNClu      = 3;

% Calculates the centroid of ward clusters
in1CluWard  = cluster(db2Tree, 'maxclust', inNClu);
in1CluWard  = SortClusterBySize(in1CluWard); %Sort the clusters to ensure consistency
db2CtrWard  = zeros(inNClu, size(db2ZData_GAD, 2));
for iClu = 1:inNClu
    db2CtrWard(iClu,:) = mean(db2ZData_GAD(in1CluWard == iClu, :));
end

% Does the K-Mean correction and prints the silhouette metrics
in1CluKM    = kmeans(db2ZData_GAD, inNClu, 'start', db2CtrWard);
%% Initialize figure and table variables
hFIG = []; iFig = 0; cFIG_LABEL = {};
% cTAB = {}; iTab = 0; cTAB_LABEL = {};% Print the list of the cells in each cluster

%Loops over clusters
for iClu = 1:inNClu
    % Selects the fraction of the data
    bl1Clu          = in1CluKM == iClu;
    db2Data_Clu     = db2Data_GAD(bl1Clu, cPAR_SEL{2});
    db2ZData_Clu    = zscore(db2Data_Clu);

    % Clusters cell according to Ephy or molecular markers
    % Clusters the data with ward methods
    db2Tree     = linkage(db2ZData_Clu, 'ward');
%     db2Inconst  = inconsistent(db2Tree); % Use the inconsistency to determine the number of cluster
%     db1Inconst  = db2Inconst(end:-1:1, 4); %Inverts the order of the rows
%     inNClu      = max(find(db1Inconst < 1, 1, 'first'), 2); % Find the first node with inconsistency inferior to one
    db1BIC = ClusterBIC(db2ZData_Clu, db2Tree, min(10, sum(bl1Clu)));
%     [~, inNSubClu] = min(db1AIC);
    db1MinD = ClusterDist(db2Tree, min(10, sum(bl1Clu)));
    [~, inNSubClu] = max(db1MinD);
    
    % Plots the dendrogram
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Dendrogram_SubClusterEphy_Clu' num2str(iClu)]);
    subplot(4, 1, [1 2])
    if inNSubClu == 1
        dendrogram(db2Tree, size(db2ZData_Clu, 1));
    else
        dendrogram(db2Tree, size(db2ZData_Clu, 1), ...
            'ColorThreshold', median(db2Tree(end - inNSubClu + 1:end - inNSubClu + 2, 3)));
    end
    title(sprintf('Subclusters Ephy Cluster%d', iClu));
    subplot(4, 1, 3)
    plot(1:length(db1BIC), db1BIC, 'ko--'); ylabel('BIC'), xlabel('Number of clusters'), xlim([0 11]);
    subplot(4, 1, 4)
    plot(1:length(db1BIC), db1MinD, 'ko--'); ylabel('Shortest Branch Length'), xlabel('Number of clusters'), xlim([0 11]);
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