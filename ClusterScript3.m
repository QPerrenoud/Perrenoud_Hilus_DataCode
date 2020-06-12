clear, close all
chInputFolder   = 'C:\Users\perre\Dropbox\Leclerc et al. 2018\Brain Structure Function\New Analysis';
% chInputFolder   = '/home/quentin/Dropbox/Leclerc et al. 2018/Brain Structure Function/New Analysis';
% chInputFolder   = '/Users/Quentin/Dropbox/Leclerc et al. 2018/Brain Structure Function/New Analysis';
cd(chInputFolder);
addpath(genpath('Utilities'))
chInputFile     = 'Input_Matrix_NewClustering';
chOutputFolder  = 'Figures';
chSubFolder     = 'GAD_PCALoadingTest';
% db1SizFig       = [100, 100, 1600, 900];
db1SizFig       = [100, 100, 800, 500];

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
db1PltQ = [.025; .25; .5; .75; .975];
iMed = find(db1PltQ == .5);
cQ_Color = {[.7 .7 .7], [.4 .4 .4], [.4 .4 .4], [.7 .7 .7]};
inNQ = length(db1PltQ);

% cTAB = {}; iTab = 0; cTAB_LABEL = {};% Print the list of the cells in each cluster

% Loops through condition
for iPSl = 1:inNSel

    % Calculate the significance and zscore of pca for each cluster
    warning off
    fprintf('%s:\r', cCLU_MSEL{iPSl});
    [db1PVal_Var, db1Z_Var] = deal(nan(1, inNClu));
    db2QLoadVar = nan(inNQ, inNClu);
    for iClu = 1:inNClu
        db2Data_CluPSel     = db2ZData_GAD(in1CluKM == iClu, cPAR_SEL{iPSl});
        db2ZData_CluPSel    = zscore(db2Data_CluPSel);
        [db1PVal_Var(iClu), db1Z_Var(iClu), db1LoadVar] = ...
            PCAPermTest(db2ZData_CluPSel);
        db2QLoadVar(:, iClu) = quantile(zscore(db1LoadVar), db1PltQ);
        fprintf('Cluster %d: Z_Var: %.2f (p = %.3f)\r', ....
            iClu, db1Z_Var(iClu), db1PVal_Var(iClu))
    end
    warning on
    %     pause(0.1), %close all
    
    % Plots a figure of the test results
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['LoadingVariance_' cCLU_MSEL{iPSl}]);
    for iClu = 1:inNClu
        plot(iClu, db1Z_Var(iClu), 'o', 'Color', cCOLOR{iClu}); hold on
        db1Y = [-.25 .25 .25 -.25] + iClu;
        for iQ = 1:inNQ - 1
            fill(db1Y, [[1 1] * db2QLoadVar(iQ, iClu) [1 1] * db2QLoadVar(iQ + 1, iClu)] , ...
                cQ_Color{iQ}, 'LineStyle', 'none', 'FaceAlpha', .3)
        end
        if ~isempty(iMed), plot([-.25 .25] + iClu, [1 1] * db2QLoadVar(iMed, iClu), 'r'); end
    end
    db1SigLevel = [.05 .01 .001];
    for iLvl = 1:length(db1SigLevel)
        db1SigY = db1Z_Var + (max(db1Z_Var) * .05) * iLvl;
        db1SigY(db1PVal_Var > db1SigLevel(iLvl)) = NaN;
        plot(1:inNClu, db1SigY, 'k*')
    end
    xlabel('Cluster'), ylabel('Loading Variance (Z-Score H0)')
    title(cCLU_MSEL{iPSl})
end
%% Plots an exemple of the PCA loading procedure
iPSl = 1;
iClu = 2;

% Calculate the significance and zscore of pca for the cluster 3
db2Data_CluPSel     = db2ZData_GAD(in1CluKM == iClu, cPAR_SEL{iPSl});
db2ZData_CluPSel    = zscore(db2Data_CluPSel);
warning off
[~, db1Z_Var, db1LoadVar] = ...
    PCAPermTest(db2ZData_CluPSel);
warning on

% Plots a histogram of the distrubution and an exemple threshold to
% illustrate the procedure
iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
cFIG_LABEL = cat(1, cFIG_LABEL, ['ProcExemple_' cCLU_MSEL{iPSl}]);
[db1H, db1Bin] = hist(zscore(db1LoadVar), 50);
bar(db1Bin, db1H, 'LineStyle', 'none'); hold on
db1YL = ylim;
plot([1 1] * db1Z_Var, db1YL, 'r')
ylabel('Permutation'), xlabel('Loading Variance (Z-Score H0)')
title(sprintf('Exemple procedure, %s parameters, Cluster %d', cCLU_MSEL{iPSl}, iClu))
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