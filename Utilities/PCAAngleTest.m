function [dbPC_DotP_Test, dbPVal, db1CI, db2Score_H0, db2Score_Test, ...
    db1VarExp_H0, db1VarExp_Test] = PCAAngleTest(db2Mat_H0, db2Mat_Test, inNIter)

%Deals with optional arguments
narginchk(1, 3)
if nargin < 3, inNIter = 10000; end

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
db1Mu_H0    = mean(db2Mat_H0);
db1SD_H0    = std(db2Mat_H0);
db2Z_H0     = (db2Mat_H0 - db1Mu_H0)./db1SD_H0;
db2Z_Test   = (db2Mat_Test - db1Mu_H0)./db1SD_H0;

% calculates the ref loading of the first pca;
[db2PC_Coef, ~, db1VarExp_H0] = pca(db2Z_H0);
db1VarExp_H0    = db1VarExp_H0 ./ sum(db1VarExp_H0);
db1PC_Ref       = db2PC_Coef(:, 1)'; 
db2Score_H0     = db2Z_H0 * db2PC_Coef; 
db2Score_Test   = db2Z_Test * db2PC_Coef;

%Initializes distribution variable
db2PC_Coef = nan(inNPar, inNIter);

%Computes the coordinates of the first PCA axis on random subsample of H0
%of the same size as test
for iItr = 1:inNIter
    db2Z_H0_Itr         = db2Z_H0(randperm(inNObs_H0, inNObs_Tst), :);
    db2PC_Coef_Itr      = pca(db2Z_H0_Itr);
    db2PC_Coef(:, iItr) = db2PC_Coef_Itr(:, 1);
end
db1PC_DotP_H0   = abs(db1PC_Ref * db2PC_Coef); % Dot product to the average (= to the cosine of the angle to ref since PC are unit vectors);
 
%Now calculates the test dot product and derives a p-value
[db2PC_Coef_Test, ~, db1VarExp_Test] = pca(db2Z_Test);
db1VarExp_Test  = db1VarExp_Test ./ sum(db1VarExp_Test);
dbPC_DotP_Test  = abs(db1PC_Ref * db2PC_Coef_Test(:, 1));
dbPVal          = mean(db1PC_DotP_H0 < dbPC_DotP_Test);
db1CI           = quantile(db1PC_DotP_H0, [.05 .25 .5 .75 1])';