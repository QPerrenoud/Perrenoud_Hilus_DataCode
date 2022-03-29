clear, close all
chInputFolder   = 'D:\Perrenoud_Leclerc_etal\Perrenoud_Hilus_DataCode'; % < ---- CHANGE TO LOCAL PATH
cd(chInputFolder);
addpath(genpath('Utilities'))
chInputFile     = 'Input_Matrix_Hilus';
chOutputFolder  = 'Figures';
chSubFolder     = '2__ParameterSets_Fig2f';
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
% Opens an output file
warning off
mkdir(chOutputFolder);
mkdir(fullfile(chOutputFolder, chSubFolder));
warning on
hFID = fopen(fullfile(chOutputFolder, chSubFolder, 'Output.txt'), 'w');

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
    db2Data_PSel_1      = db2Data(cSEL{iPop}, cPAR_SEL{1});
    db2ZData_PSel_1     = zscore(db2Data_PSel_1);
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
    db1Sil_1    = silhouette(db2ZData_PSel_1, in1CluKM);
    db1Sil_2    = silhouette(db2ZData_PSel, in1CluKM);
    fprintf(hFID, '%s : \r', cCLU_MSEL{iPSl});
    fprintf(hFID, '*************************\r\r');
    fprintf(hFID, '-------- Silhouette -----------------------------\r');
    fprintf(hFID, '--Original parameter set\r');
    fprintf(hFID, '%d cell out %d have a negative silhouette value\r', sum(db1Sil_1 < 0), inNObs);
    fprintf(hFID, 'Mean silhouette: %.3f\r', mean(db1Sil_1));
    fprintf(hFID, 'Restricted parameter set\r');
    fprintf(hFID, '%d cell out %d have a negative silhouette value\r', sum(db1Sil_2 < 0), inNObs);
    fprintf(hFID, 'Mean silhouette: %.3f\r\r', mean(db1Sil_2));
    
    %Plots the silhouette
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Silhouette_' cCLU_MSEL{iPSl}]);
    silhouette(db2ZData_PSel, in1CluKM);
    
    % Plost the value of electrophysiological parameters
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['ParEphy_' cCLU_MSEL{iPSl}]);
    fprintf(hFID, '-------- Ephys ---------------------------------\r');
    for iPar = in1Phy
        fprintf(hFID, '%s: ++++++++++++++++++++++++++++++++\r', cPAR{iPar});
        subplot(4, 4, iPar), hold on
        for iClu = 1:inNClu
            bl1CluCells = in1CluKM == iClu;
            inNCells    = sum(bl1CluCells);
            dbMean  = mean(db2Data(bl1CluCells, iPar));
            dbSEM   = sqrt(var(db2Data(bl1CluCells, iPar))/inNCells);
            errorbar(iClu, dbMean, dbSEM, dbSEM, 'o', 'Color', cCOLOR{iClu}, 'LineWidth', 2);
            fprintf(hFID, 'Cluster %d: %.3f (+- %.3f)\r', iClu, dbMean, dbSEM); % Print results to file
        end
        set(gca, 'XTick', 1:inNClu);
        xlim([0 inNClu + 1]);
        xlabel('Clusters');
        title(cPAR(iPar));
        
        % Prints the significance to test 
        fprintf(hFID, '--------------Sigficance (Mann Witney):\r');
        for iClu = 1:inNClu, fprintf(hFID, '\tC%d', iClu); end
        for iCl1 = 1:inNClu - 1
            fprintf(hFID, '\r');
            fprintf(hFID, 'C%d', iCl1); for iClu = 1:iCl1, fprintf(hFID, '\t'); end
            for iCl2 = iCl1 + 1:inNClu
                dbPVal = ranksum(db2Data(in1CluKM == iCl1, iPar), db2Data(in1CluKM == iCl2, iPar));
                fprintf(hFID, '\t%.3f', dbPVal);
            end
        end
        fprintf(hFID, '\r');
        fprintf(hFID, '\r');
    end
    
    % Plots the marker histograms of the clusters
    iFig = iFig + 1; hFIG(iFig) = figure('Position', db1SizFig);
    cFIG_LABEL = cat(1, cFIG_LABEL, ['Histogram_' cCLU_MSEL{iPSl}]);
    cMARKER     = {'VGluT1', 'GAD65', 'GAD67', 'NOS1', 'PV', 'CR', 'NPY', 'SOM', 'CCK'};
    in1MrkIdx   = cellfun(@(x) find(ismember(cLABELS, x)), cMARKER);
    db2DatMark  = db2DataAll_PSel(:, in1MrkIdx);
    db2DatMark(:, 2) = db2DatMark(:, 2) | db2DatMark(:, 3);
    cMARKER(3) = []; cMARKER{2} = 'GAD';
    db2DatMark(:, 3) = [];
    fprintf(hFID, '-------- Molecular markers -----------------------------\r');
    fprintf(hFID, '\t');
    for iMrk = 1:length(cMARKER), fprintf(hFID, '%s\t', cMARKER{iMrk}); end
    fprintf(hFID, '\r');
    if inNClu <= inNPltCol, inRow = 1; inCol = inNClu;
    else, inRow = ceil(inNClu/inNPltCol); inCol = inNPltCol; end
    [db2Prop, in2NObs] = deal(nan(inNClu, length(cMARKER))); 
    for iClu = 1:inNClu
        subplot(inRow, inCol, iClu)
        db2Prop(iClu, :)    = mean(db2DatMark(in1CluKM == iClu, :), 1);
        in2NObs(iClu, :)     = sum(in1CluKM == iClu); 
        bar(1:length(cMARKER), db2Prop(iClu, :))
        set(gca, 'XTick', 1:length(cMARKER), 'XTickLabel', cMARKER')
        ylim([0 1.2])
        title(sprintf('Cluster %d: N = %d', iClu, sum(in1CluKM == iClu)));
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
    for iCl1 = 1:inNClu - 1
        for iCl2 = iCl1 + 1:inNClu
            db1PVal = BinomialTest2(in2NObs(iCl1, :), db2Prop(iCl1, :), in2NObs(iCl2, :), db2Prop(iCl2, :));
            fprintf(hFID, 'C%d-C%d', iCl1, iCl2);
            for iMrk = 1:length(cMARKER), fprintf(hFID, '\t%.3f', db1PVal(iMrk)); end
            fprintf(hFID, '\r');
        end
    end
    fprintf(hFID, '\r');
    fprintf(hFID, '\r');
    
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

% Print number of reassigned cells
fprintf(hFID, '\r');
fprintf(hFID, '*******************************************************************\r');
fprintf(hFID, '-------- Comparison of parameter sets -----------------------------\r');
for iSel = 2:length(cCLU_MSEL)
    inNDif = sum(in2_CluKMean(:, 1) ~= in2_CluKMean(:, iSel));
    fprintf(hFID, '%s:\t%d reassignment\r', cCLU_MSEL{iSel}, inNDif);
end

% Close the output writing file
fclose(hFID);
%%
save(fullfile(chOutputFolder, chSubFolder, 'Workspace'), '-v7.3')
savefig(hFIG,fullfile(chOutputFolder, chSubFolder, 'Figures'));

SaveFig(fullfile(chOutputFolder, chSubFolder), hFIG, cFIG_LABEL);
close all