function [db1Dist] = ClusterDist(db2Tree, inNCluMax)

% Determines the number of nodes and observation
inNNode = size(db2Tree, 1);
inNObs  = inNNode + 1;

db1Dist = zeros(1, min(inNNode, inNCluMax - 1));
for iNde = 1:length(db1Dist - 1)
    % We start a the topmost node which is the last of the tree and compute
    % the height of the node
    inIdx       = inNNode - iNde + 1; 
    dbD_Nde     = db2Tree(inIdx, 3);
    
    % We then find the height of the subnodes at the end each branch of the node
    inIdx_B1    = db2Tree(inIdx, 1); % Index of the node at the end of branch 1;
    if inIdx_B1 > inNObs, dbD_B1 = db2Tree(inIdx_B1 - inNObs, 3); 
    else, dbD_B1 = 0; end % Height of the node
    inIdx_B2    = db2Tree(inIdx, 2); % Index of the node at the end of branch 2;
    if inIdx_B2 > inNObs, dbD_B2 = db2Tree(inIdx_B2 - inNObs, 3); 
    else, dbD_B2 = 0; end % Height of the node
    
%     % The overall distance of the node is the the sum of the distance of
%     % each branch
%     db1Dist(iNde) = 2 * dbD_Nde - dbD_B1 - dbD_B2;

    % Alternatively takes the minimum branch length (Thorndike procedure)
    db1Dist(iNde + 1) = min(dbD_Nde - dbD_B1, dbD_Nde - dbD_B2);
end