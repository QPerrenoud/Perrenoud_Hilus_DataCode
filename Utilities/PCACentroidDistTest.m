function [dbPVal] = PCACentroidDistTest(db2Mat_H0, db2Mat_Test, inNIter, inNPC)

%Deals with optional arguments
narginchk(1, 3)
if nargin < 3, inNIter = 10000; end
if nargin < 4, inNPC = 2; end

% Checks that the input is properly formatted
[inNObs_H0, inNPar] = size(db2Mat_H0);
[inNObs_Tst, inNPar_Tst] = size(db2Mat_Test);
if inNPar_Tst ~= inNPar
    error('db2MatH0 and db2MatTest must have the same number of columns');
end
if inNObs_Tst > inNObs_H0/2
    error('db2MatH0 must have at least twice as many rows (observations) as db2MatTest');
end
if mod(inNIter, 1) ~= 0 || inNIter <= 0
    error('inNIter must be a positive integer')
end

% Initialize the random number generator for reproducibility
rng(1949);

% zscore parameters to their distribution accross H0 and test samples
db2CatMat   = cat(1, db2Mat_H0, db2Mat_Test);
db2ZScore   = zscore(db2CatMat);
db2Z_H0     = db2ZScore(1:inNObs_H0, :);
db2Z_Test   = db2ZScore(inNObs_H0 + 1: end, :);

% calculates the ref loading of the first pca;
db2PC_Coef      = pca(db2Z_H0);
db2Score_H0     = db2Z_H0 * db2PC_Coef;
db2Score_Test   = db2Z_Test * db2PC_Coef;

% calculates the ref loading of the first pca;
db1Cntrd_H0  = mean(db2Score_H0(:, 1:inNPC), 1);

%Initializes distribution variable
db2Cntrd_Itr  = nan(inNIter, inNPC);

%Computes the coordinates of the first PCA axis on random subsample of H0
%of the same size as test
for iItr = 1:inNIter
    db2Score_H0_Itr         = db2Score_H0(randperm(inNObs_H0, inNObs_Tst), :);
    db2Cntrd_Itr(iItr, :)   = mean(db2Score_H0_Itr(:, 1:inNPC), 1);
end
%db1Cntrd_H0  = mean(db2Cntrd_Itr, 1);
db1Dist_H0   = sqrt(sum((db2Cntrd_Itr - db1Cntrd_H0).^2, 2));
 
%Now calculates the test dot product and derives a p-value
db1Cntr_Test    = mean(db2Score_Test(:, 1:inNPC), 1);
dbDist_Test     = sqrt(sum((db1Cntr_Test - db1Cntrd_H0).^2, 2));
dbPVal          = mean(db1Dist_H0 > dbDist_Test);