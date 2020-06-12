function [dbPVal, dbZVar, db1LoadVar] = PCAPermTest(db2Mat, inNIter)
%DBPVAL = PCAPermTest(DB2MAT, [INNITER], [INNCMP])
% %Estimates the p-value DBPVAL of the variance explained by the first
% INNCMP components of a PCA performed on the matrice DB2MAT where rows
% are observation and columns are parameters. The p-value is computed using
% a permutation %test where the order of the value of each parameters is
% scrambled along the columns of DB2MAT. INNITER controls the number of
% permutations used to estimate the distribution.

%Deals with optional arguments
narginchk(1, 2)
if nargin < 2, inNIter = 10000; end

%Estimates the size of the matrix
[inNObs, inNVar] = size(db2Mat);

%Loops through iterations
db1LoadVar = deal(nan(1, inNIter));
for iItr = 1:inNIter
    db2MatPerm = db2Mat;
    for iVar = 1:inNVar
        db2MatPerm(:, iVar) = db2MatPerm(randperm(inNObs), iVar);
    end
    [~, ~, ~, ~, db1PCAVar] = pca(db2MatPerm);
    db1LoadVar(iItr) = var(db1PCAVar);
end

%Estimates the p-value
[~, ~, ~, ~, db1PCAVar] = pca(db2Mat);
dbVarVar    = var(db1PCAVar);
dbPVal      = mean(db1LoadVar > dbVarVar);
dbZVar      = (dbVarVar - mean(db1LoadVar))./std(db1LoadVar);