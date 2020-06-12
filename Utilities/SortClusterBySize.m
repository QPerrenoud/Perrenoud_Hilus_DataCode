function in1CluIdxSorted = SortClusterBySize(in1CluIdxRaw)
%Function to reassign cluster number as a function of size. This can ensure
%consistency in cluster display

% Gets the number of cluster
inNClu = max(in1CluIdxRaw);

% Initializes the sorted cluster indices
in1CluIdxSorted = zeros(size(in1CluIdxRaw));

% Estimates the number of observation in each cluster
db1CluNObs      = hist(in1CluIdxRaw, 1:inNClu);
[~, in1SrtIdx]   = sort(db1CluNObs);
for iClu = 1:inNClu
    in1Idx = in1CluIdxRaw == in1SrtIdx(iClu);
    in1CluIdxSorted(in1Idx) = iClu;
end