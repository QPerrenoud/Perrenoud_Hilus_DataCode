clear, close all
chInputFolder   = 'C:\Users\perre\Dropbox\Leclerc et al. 2018\Brain Structure Function\New Analysis';
% chInputFolder   = '/home/quentin/Dropbox/Leclerc et al. 2018/Brain Structure Function/New Analysis';
% chInputFolder   = '/Users/Quentin/Dropbox/Leclerc et al. 2018/Brain Structure Function/New Analysis';
cd(chInputFolder);
addpath(genpath('Utilities'))
chInputFile     = 'Input_Matrix_NewClustering';
chOutputFolder  = 'Figures';
chSubFolder     = 'Global_vs_GAD';
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
db2ZData    = zscore(db2Data);

% Seed random number generator
rng(1984)
%% Clusters GABAergic cells based on electrophysiological or molecular marker
% Defines sub populations of interest for the analysis
cCLU_POP    = {'Global', 'GAD'};
cSEL        = {true(size(db2ZData(:, 1))) db2Data(:, 18) == 1};

% Initializes a color vector to plot clusters in different colors
cCOLOR = {[255 222 23]./255, [57 181 74]./255, [185 82 159]./255, [.4 .8 0], [0 .5 0], [.9 .45 0], [.9 0 0], ...
    [.9 0 .45], [0 0 .9], [0 .45 .9], [.45 0 .9], [.45 .9 0], [0 .9 .45]};
%% CLusters GABAergic neurons
% Selects the fraction of the data
iPop = 2;
db2Data_GAD         = db2Data(cSEL{iPop}, :);
db2ZData_GAD        = zscore(db2Data_GAD);
cNEU_GAD            = cNEURONS(cSEL{iPop});

% Clusters cell according to Ephy or molecular markers
% Clusters the data with ward methods
db2Tree_GAD     = linkage(db2ZData_GAD, 'ward');
inNClu_GAD      = 3;

% Calculates the centroid of ward clusters
in1CluWard_GAD  = cluster(db2Tree_GAD, 'maxclust', inNClu_GAD);
in1CluWard_GAD  = SortClusterBySize(in1CluWard_GAD); %Sort the clusters to ensure consistency
db2CtrWard_GAD  = zeros(inNClu_GAD, size(db2ZData_GAD, 2));
for iClu = 1:inNClu_GAD
    db2CtrWard_GAD(iClu,:) = mean(db2ZData_GAD(in1CluWard_GAD == iClu, :));
end

% Does the K-Mean correction and prints the silhouette metrics
in1CluKM_GAD    = kmeans(db2ZData_GAD, inNClu_GAD, 'start', db2CtrWard_GAD);
%% Print the list of the cells in each cluster
% Find the cluster of CL427 (Mossy neurons, should not be here); CL430 and CL477
% in1427 = find(ismember(cNEU_SELL, 'CL427'));
% inClu427 = in1CluKM(in1427);
% in1430 = find(ismember(cNEU_SELL, 'CL430'));
% inClu430 = in1CluKM(in1430);
% in1477 = find(ismember(cNEU_SELL, 'CL477'));
% inClu477 = in1CluKM(in1477);

% Print the list of the cells in each cluster
for iClu = 1:inNClu_GAD
    in1CluObs = find(in1CluKM_GAD == iClu);
    cOBS = cNEU_GAD(in1CluObs);
    db2DataSel = db2Data_GAD(in1CluObs, :);
    fprintf('Cluster %d:\r', iClu)
    fprintf('Cell:\t\t')
    for iPar = 19:24, fprintf('%s\t\t', cPAR{iPar}(1:min(end, 3))); end
    fprintf('\r')
    for iObs = 1:length(cOBS)
        fprintf('%s\t\t', cOBS{iObs});
        for iPar = 19:24, fprintf('%d\t\t', db2DataSel(iObs, iPar)); end
        fprintf('\r')
    end
    fprintf('\r');
end

% Notes for potentials exemple cells
%Cluster 1: 
%   CL209, CL221, CL260, CL300                      -> CL300
%
%Cluster 2:
%   - NOS1: CL43 CL290 CL307                        -> CL307
%   - PV:   CL43 CL73 CL179 CL449                   -> CL73
%   - CR:   CL55 CL99 CL104 CL105 CL179 CL245       -> CL245
%   - NPY:  CL178 CL313                             -> CL313
%% CLusters All neurons, plot the dendrogram and attempt to visualize where GABAergic neurons of different cluster fall
% Initialize the figure
clear hFIG
iFig = 0; cFIG_LABEL = {};

% Selects the fraction of the data
iPop = 1;
db2Data_All         = db2Data(cSEL{iPop}, :);
db2ZData_All        = zscore(db2Data_All);
cNEU_ALL            = cNEURONS(cSEL{iPop});

% Clusters cell according to Ephy or molecular markers
% Clusters the data with ward methods
db2Tree_All     = linkage(db2ZData_All, 'ward');
inNClu_All      = 2;

% Plots the dendrogram
iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
cFIG_LABEL = cat(1, cFIG_LABEL, ['Dendrogram_All_GAD_Cells']);
[~, ~, in1WrdIdx_All] = dendrogram(db2Tree_All, size(db2ZData_All, 1),....
    'ColorThreshold', median(db2Tree_All(end - inNClu_All + 1:end - inNClu_All + 2, 3)));
% Marks the branches corresponding to the different clusters of GABAergic
% cells idendified in the subset of the selections
hold on
for iClu = 1:inNClu_GAD
    in1CluObs = find(in1CluKM_GAD == iClu);
    in1Idx1 = find(ismember(cNEU_ALL, cNEU_GAD(in1CluObs)));
    in1Idx2 = find(ismember(in1WrdIdx_All, in1Idx1));
    plot(in1Idx2, zeros(size(in1Idx2)), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cCOLOR{iClu});
end
%% Save the figures and the tables
warning off
mkdir(chOutputFolder);
mkdir(fullfile(chOutputFolder, chSubFolder));
warning on)

for iFig = 1:length(hFIG)
    exportfig(hFIG(iFig) ,fullfile(chOutputFolder, chSubFolder, cFIG_LABEL{iFig}), 'Color', 'rgb', 'Renderer', 'painters');
    set(hFIG(iFig), 'Units', 'centimeters')
    db1Pos = get(hFIG(iFig), 'Position');
    set(hFIG(iFig), 'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', 'PaperSize', [db1Pos(3) db1Pos(4)])
    print(hFIG(iFig), fullfile(chOutputFolder, chSubFolder, cFIG_LABEL{iFig}), '-dpdf')
end