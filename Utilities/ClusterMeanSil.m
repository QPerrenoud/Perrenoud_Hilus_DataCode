function [db1MeanSil] = ClusterMeanSil(db2Data, db2Tree, inNCluMax)

% Initializes the output variable
db1MeanSil  = zeros(1, inNCluMax);

for iNCl = 2:inNCluMax
    % Calculates the centroid of ward clusters
    in1CluWard  = cluster(db2Tree, 'maxclust', iNCl);
    db2CtrWard  = zeros(iNCl, size(db2Data, 2));
    for iClu = 1:iNCl
        db2CtrWard(iClu,:) = mean(db2Data(in1CluWard == iClu, :));
    end
    
    % Does the K-Mean correction
    in1CluKM  = kmeans(db2Data, iNCl, 'start', db2CtrWard , 'onlinephase', 'off');
    
    % Calculates the Bayesian information criterion of the clustering
    db1MeanSil(iNCl) = mean(silhouette(db2Data, in1CluKM));
end