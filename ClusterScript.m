clear, close all
chInputFolder   = 'C:\Users\perre\Dropbox\Leclerc et al. 2018\Brain Structure Function\New Analysis';
% chInputFolder   = '/home/quentin/Dropbox/Leclerc et al. 2018/Brain Structure Function/New Analysis';
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

% Construct the parameter matrix
in1Idx      = [1:4 7 9 13:16 19:20 24 28 39 40 41:42 44 46:48 50:51]; %Loads the parameters
cPAR        = cLABELS(in1Idx); %
db2Data     = db2DataAll(:, in1Idx);
cPAR{18}    = 'GAD';
db2Data(:, 18) = db2DataAll(:,42) | db2DataAll(:, 43); % Merges the expression of the two gads
db2ZData    = zscore(db2Data);

% Seed random number generator
rng(1984)
%% Plots clusters and parameters
% Defines sub populations of interest for the analysis
cCLU_POP    = {'Global', 'GAD'};
cSEL        = {true(size(db2ZData(:, 1))) db2Data(:, 18) == 1};

% Initialize figure and table variables
hFIG = []; iFig = 0; cFIG_LABEL = {}; inNPltCol = 3;
cTAB = {}; iTab = 0; cTAB_LABEL = {};

% Loops through populations of interest
for iPop = 1:length(cCLU_POP)
    % Selects the fraction of the data
    %db2ZData_Pop    = db2ZData(cSEL{iPop}, :); %Out 2020-02-23
    db2Data_Pop     = db2Data(cSEL{iPop}, :); % Added 2020-02-23
    db2ZData_Pop    = zscore(db2Data_Pop);
    db2DataAll_Pop  = db2DataAll(cSEL{iPop}, :);
    [inNObs, inNVar] = size(db2ZData_Pop);
    
    
    % Clusters the data with ward methods
    db2Tree     = linkage(db2ZData_Pop, 'ward');
%     db2Inconst  = inconsistent(db2Tree); % Use the inconsistency to determine the number of cluster
%     db1Inconst  = db2Inconst(end:-1:1, 4); %Inverts the order of the rows
%     [~, inNClu] = min(db1Inconst(1:10));
%     inNClu      = find(db1Inconst < 1, 1, 'first'); % Find the first node with inconsistency inferior to one
    db1BIC      = ClusterBIC(db2ZData_Pop, db2Tree, 10);
%     db1AIC      = ClusterLogLikelihood(db2ZData_Pop, 10);
%     [~, inNClu] = min(db1AIC);
    db1MinD = ClusterDist(db2Tree, 10);
    [~, inNClu] = max(db1MinD);

    % Plots the dendrogram
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Dendrogram_' cCLU_POP{iPop}]);
    subplot(4, 1, [1 2])
    dendrogram(db2Tree, size(db2ZData_Pop, 1), 'ColorThreshold', median(db2Tree(end - inNClu + 1:end - inNClu + 2, 3)));
    subplot(4, 1, 3)
    plot(1:10, db1BIC, 'ko--'); ylabel('BIC'), xlabel('Number of clusters'), xlim([0 11]);
    subplot(4, 1, 4)
    plot(1:10, db1MinD, 'ko--'); ylabel('Shortest Branch Length'), xlabel('Number of clusters'), xlim([0 11]);
    
    % Calculates the centroid of ward clusters
    in1CluWard  = cluster(db2Tree, 'maxclust', inNClu);
    % Ensures that the clusters are ordered by number of cells
    in1CluWard = SortClusterBySize(in1CluWard); 
   
    % Initializes centroids
    db2CtrWard=zeros(inNClu, size(db2ZData_Pop, 2));
    for iClu = 1:inNClu
        db2CtrWard(iClu,:) = mean(db2ZData_Pop(in1CluWard == iClu, :));
    end
    
    % Does the K-Mean correction and prints the silhouette metrics
    in1CluKM    = kmeans(db2ZData_Pop, inNClu, 'start', db2CtrWard);
    db1Sil      = silhouette(db2ZData_Pop, in1CluKM);
    fprintf('%d cellules sur %d ont une valeur de silhouette negative\r', length(db1Sil(db1Sil < 0)), inNObs);
    fprintf('Silhouette moyenne: %.3f\r', mean(db1Sil));
    % Creates the table of correspondance between k-mean and ward
    db1Tab = nan(inNClu + 1);
    for iClW = 1:inNClu
        for iClK = 1:inNClu
            db1Tab(iClW, iClK) = sum(in1CluWard == iClW & in1CluKM == iClK);
        end
    end
    db1Tab(:, end) = nansum(db1Tab, 2);
    db1Tab(end, :) = nansum(db1Tab, 1);
    iTab = iTab + 1;
    cTAB_LABEL = cat(1, cTAB_LABEL, sprintf('Ward_vs_KMeans_%s', cCLU_POP{iPop}));
    Clusters = [cellfun(@(x) sprintf('Clu%d', x), num2cell(1:inNClu), 'UniformOutput', false) 'KMeans']';
    chTabString = 'cTAB{iTab} = table(Clusters, ';
    for iClu = 1:inNClu
        eval(sprintf('Clu%d = db1Tab(:, %d);', iClu, iClu));
        chTabString = cat(2, chTabString, sprintf('Clu%d, ', iClu));
    end
    Ward = db1Tab(:, end);
    chTabString = cat(2, chTabString, 'Ward);');
    eval(chTabString)
    disp(cTAB{iTab})

    %Plots the silhouette
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Silhouette_' cCLU_POP{iPop}]);
    silhouette(db2ZData_Pop, in1CluKM);
        
    % Plots the marker histograms of the clusters
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Histogram_Mol_' cCLU_POP{iPop}]);
    cMARKER     = {'VGluT1', 'GAD65', 'GAD67', 'NOS1', 'CA', 'PV', 'CR', 'NPY', 'VIP', 'SOM', 'CCK'};
    in1MrkIdx   = cellfun(@(x) find(ismember(cLABELS, x)), cMARKER);
    db2DatMark  = db2DataAll_Pop(:, in1MrkIdx);
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
    pause(.1)
    
    % Plots the clusters on the two principal components
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['PCA_' cCLU_POP{iPop}]);
    [~, db1PC] = pca(db2ZData_Pop);
    gscatter(db1PC(:,1), db1PC(:,2), in1CluKM);
    
    % Scrambles variables one by one to see their effect on the mean
    % silhouette value
    inNIter = 1000; % number of iterations
    db2Scrbl = zeros(inNIter, inNVar);
    for iPar = 1:inNVar
        for iItr = 1:inNIter
            db2DatRnd           = db2ZData_Pop;
            db2DatRnd(:,iPar)   = db2DatRnd(randperm(inNObs), iPar);
            in1CluKM_Rnd        = kmeans(db2DatRnd, inNClu, 'start', db2CtrWard);
            db2Scrbl(iItr,iPar) = mean(silhouette(db2DatRnd, in1CluKM_Rnd));
        end
    end
    db1SigThr   = quantile(db2Scrbl, .95);
    bl1Sig      = db1SigThr < mean(db1Sil);
    
    %Plots the result of the scrambling
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Scrambling_' cCLU_POP{iPop}]);
    bar(1:inNVar, mean(db2Scrbl)), hold on
    errorbar(1:inNVar, mean(db2Scrbl), nan(1, inNVar), std(db2Scrbl), '.k');
    ylabel('Silhouette (mean + s.d.)');
    set(gca, 'XTick', 1:inNVar, 'XTickLabel', cPAR);
    xtickangle(60)
    title('Effects of the scrampbling on clustering quality')
    db1SigPlt = nan(1, inNVar); db1SigPlt(bl1Sig) = 1;
    plot(1:inNVar, db1SigPlt * 1.1 * mean(db1Sil), '*k')
    hLINE = plot(xlim, [1 1]*mean(db1Sil), '--r');
    legend(hLINE, {'Original Silhouette'})
    ylim([0 mean(db1Sil) * 1.2]);
    
    
    % Test the association of marker within clusters
    % FrequencyTest(db2DataAll(in1Clu == inU_Clu(1), in1MrkIdx), db2DataAll(in1Clu == inU_Clu(2), in1MrkIdx));
    for iClu = 1:inNClu
        db2_PVal    = AssociationTest(db2DatMark(in1CluKM == iClu, :));
        db2_Sig     = double(db2_PVal < .05);
        db2_PVal(diag(true(1, size(db2_Sig, 2)))) = NaN;
        db2_Sig(diag(true(1, size(db2_Sig, 2)))) = NaN;
        
        % Generates and prints a table of marker association
        iTab = iTab + 1;
        cTAB_LABEL = cat(1, cTAB_LABEL, sprintf('MarkerAssociation_Pval_%s_Cluster_%d', cCLU_POP{iPop}, iClu)); 
        chTabString = 'cTAB{iTab} = table(cMARKER'', ';
        for iPar = 1:length(cMARKER)
            eval(sprintf('%s = db2_PVal(:, %d);', cMARKER{iPar}, iPar));
            if iPar < length(cMARKER), chTabString = cat(2, chTabString, sprintf('%s, ', cMARKER{iPar}));
            else, chTabString = cat(2, chTabString, sprintf('%s);', cMARKER{iPar})); end
        end
        eval(chTabString);
        fprintf('Cluster %d\r', iClu)
        disp(cTAB{iTab})
        
                % Generates and prints a table of marker association
        iTab = iTab + 1;
        cTAB_LABEL = cat(1, cTAB_LABEL, sprintf('MarkerAssociation_%s_Cluster_%d', cCLU_POP{iPop}, iClu)); 
        chTabString = 'cTAB{iTab} = table(cMARKER'', ';
        for iPar = 1:length(cMARKER)
            eval(sprintf('%s = db2_Sig(:, %d);', cMARKER{iPar}, iPar));
            if iPar < length(cMARKER), chTabString = cat(2, chTabString, sprintf('%s, ', cMARKER{iPar}));
            else, chTabString = cat(2, chTabString, sprintf('%s);', cMARKER{iPar})); end
        end
        eval(chTabString);
        fprintf('Cluster %d\r', iClu)
        disp(cTAB{iTab})
    end
end
%% Save the figures and the tables
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

for iTab = 1:length(cTAB)
    writetable(cTAB{iTab}, fullfile(chOutputFolder, chSubFolder, [cTAB_LABEL{iTab} '.xls']))
end