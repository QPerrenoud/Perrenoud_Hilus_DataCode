function [db1MeanSil_Nrm] = ClusterMeanSilNorm(db2Data, inNCluMax)
% Utility returning the normalized log likelihood of a clustering partition.
% The partition is compared to the average log likelihood of a partition on
% scrambled data

% Controls for optional arguments
narginchk(2, 3);
if nargin < 3, inNIter = 1000; end

% Gets the size of the data
[inNObs, inNPar] = size(db2Data);

% Initializes the output variable
db1MeanSil_Nrm  = zeros(1, inNCluMax);

% Computes a dendrogram on the data with Ward's method
db2Tree     = linkage(db2Data, 'ward');

for iNCl = 2:inNCluMax
    db1MS_H0 = nan(1, inNIter);
    for iItr = 1:inNIter
        % Randomized the data and computes a dendrogram
        db2DataRnd =  db2Data;
        for iPar = 1:inNPar
            db2DataRnd(:, iPar) = db2DataRnd(randperm(inNObs), iPar);
        end
        db2TreeRnd = linkage(db2DataRnd, 'ward');
        
        % Calculates the mean silhouette of the randomized data
        db1MS_H0(iItr) = MeanSil(db2DataRnd, db2TreeRnd, iNCl);
    end
    
    % Computes the observation probabilities
    dbMS = MeanSil(db2Data, db2Tree, iNCl);
    
    % Calculates the Bayesian information criterion of the clustering
    db1MeanSil_Nrm(iNCl) = (dbMS - mean(db1MS_H0))./std(db1MS_H0); 
end

function dbMeanSil = MeanSil(db2Data, db2Tree, iNCl)
    
% Calculates the centroid of ward clusters
in1CluWard  = cluster(db2Tree, 'maxclust', iNCl);
db2CtrWard  = zeros(iNCl, size(db2Data, 2));
for iClu = 1:iNCl
    db2CtrWard(iClu,:) = mean(db2Data(in1CluWard == iClu, :));
end
    
% Does the K-Mean correction
in1CluKM  = kmeans(db2Data, iNCl, 'start', db2CtrWard , 'onlinephase', 'off');

% Calculates the Bayesian information criterion of the clustering
dbMeanSil = mean(silhouette(db2Data, in1CluKM));