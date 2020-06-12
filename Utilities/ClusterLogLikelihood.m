function [db1NrmLgL] = ClusterLogLikelihood(db2Data, inNCluMax, inNIter)
% Utility returning the normalized log likelihood of a clustering partition.
% The partition is compared to the average log likelihood of a partition on
% scrambled data

% Controls for optional arguments
narginchk(2, 3);
if nargin < 3, inNIter = 1000; end

% Gets the size of the data
[inNObs, inNPar] = size(db2Data);

% Initializes the output variable
db1NrmLgL  = inf(1, inNCluMax);

% Computes a dendrogram on the data with Ward's method
db2Tree     = linkage(db2Data, 'ward');

for iNCl = 2:inNCluMax
    db1LgLH0 = nan(1, inNIter);
    for iItr = 1:inNIter
        % Randomized the data and computes a dendrogram
        db2DataRnd =  db2Data;
        for iPar = 1:inNPar
            db2DataRnd(:, iPar) = db2DataRnd(randperm(inNObs), iPar);
        end
        db2TreeRnd = linkage(db2DataRnd, 'ward');
        
        % Calculates the observation probabilities on the randomized data
        db1PRnd = ObsProb(db2DataRnd, db2TreeRnd, iNCl);
        
        % Computes the Log Likelihood of the randomized model
        db1LgLH0(iItr) = sum(log(db1PRnd));
    end
    
    % Computes the observation probabilities
    db1Prob = ObsProb(db2Data, db2Tree, iNCl);
    
    % Calculates the Bayesian information criterion of the clustering
    db1NrmLgL(iNCl) = (mean(db1LgLH0) - sum(log(db1Prob)))./std(db1LgLH0); 
end

function db1Prob = ObsProb(db2Data, db2Tree, iNCl)
    
% Calculates the centroid of ward clusters
in1CluWard  = cluster(db2Tree, 'maxclust', iNCl);
db2CtrWard  = zeros(iNCl, size(db2Data, 2));
for iClu = 1:iNCl
    db2CtrWard(iClu,:) = mean(db2Data(in1CluWard == iClu, :));
end

% Does the K-Mean correction
[~,  ~, ~, db2Dist]  = kmeans(db2Data, iNCl, 'start', db2CtrWard , 'onlinephase', 'off');
db2Dist = sqrt(db2Dist);

% Estimate the probability that each point is within its cluster as:
% pi = 1 - ((xi - ci)/(xi - cj)) where xi is the coordinate of
% observation i, ci is the centroid of the cluster to which i is
% assigned and cj is the centroid of the closest cluster.
db2SrtDist  = sort(db2Dist, 2);
db1Prob     = 1 - (db2SrtDist(:, 1)./db2SrtDist(:, 2));