%% Simulation for the BIC.
inNItr      = 100;
inNCluMax   = 100;

% Initializes the output variable
[db1BIC, db1AIC, db1LgL]  = deal(nan(inNItr, inNCluMax - 1));

for iItr = 1:inNItr
    db2Data = randn(150, 30); % Simulate some unicluster data
    
    % Gets the number of observation and the number of varialbes
    inNObs  = size(db2Data, 1);
    
    % Computes a dendrogram with Ward's method
    db2Tree     = linkage(db2Data, 'ward');
    
    for iNCl = 2:inNCluMax
        % Calculates the centroid of ward clusters
        in1CluWard  = cluster(db2Tree, 'maxclust', iNCl);
        db2CtrWard  = zeros(iNCl, size(db2Data, 2));
        for iClu = 1:iNCl
            db2CtrWard(iClu,:) = mean(db2Data(in1CluWard == iClu, :));
        end
        
        % Does the K-Mean correction
        [~,  ~, ~, db2Dist]  = kmeans(db2Data, iNCl, 'start', db2CtrWard , 'onlinephase', 'off');
        %db2Dist = sqrt(db2Dist);
        
        % Estimate the probability that each point is within its cluster as:
        % pi = 1 - ((xi - ci)/(xi - cj)) where xi is the coordinate of
        % observation i, ci is the centroid of the cluster to which i is
        % assigned and cj is the centroid of the closest cluster.
        db2SrtDist  = sort(db2Dist, 2);
        db1Prob     = 1 - (db2SrtDist(:, 1)./db2SrtDist(:, 2));
        
        % Calculates the Bayesian information criterion of the clustering
        db1BIC(iItr, iNCl - 1) = (iNCl * log(inNObs)) - (2 * sum(log(db1Prob)));
        db1AIC(iItr, iNCl - 1) = (iNCl * 2) - (2 * sum(log(db1Prob)));
        db1LgL(iItr, iNCl - 1) = sum(log(db1Prob));
    end
end

%% plots some figures

figure('Position', [100 100 1800 900]);
subplot(1, 3, 1); plot(mean(db1BIC), '--o'); title('mean BIC');
subplot(1, 3, 2); plot(mean(db1AIC), '--o'); title('mean AIC');
subplot(1, 3, 3); plot(mean(db1LgL), '--o'); title('mean Log Likelihoood');