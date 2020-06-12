function [db1Proj_H0, db1Proj_Test] = CCAxeProjection(db2Mat_H0, db2Mat_Test)
% Utilities returning the coordinates of a set of observation contained in
% DB2MAT_H0 and DB2MAT_TEST when projected on the axis linking the centroid
% of each data set. The coordinates are normalized so that the centroid of
% H0 and Test occupy positon 0 and 1 respectively.

%Deals with optional arguments
narginchk(1, 2)

% Checks that the input is properly formatted
[~ , inNPar] = size(db2Mat_H0);
[~, inNPar_Tst] = size(db2Mat_Test);
if inNPar_Tst ~= inNPar
    error('db2MatH0 and db2MatTest must have the same number of columns');
end

% zscore parameters to their distribution accross H0 and test samples
db1Mu_H0    = mean(db2Mat_H0);
db1SD_H0    = std(db2Mat_H0);
db2Z_H0     = (db2Mat_H0 - db1Mu_H0)./db1SD_H0;
db2Z_Test   = (db2Mat_Test - db1Mu_H0)./db1SD_H0;

% Calculates the cooerdinates of the centroids
db1Ctr_H0   = mean(db2Z_H0);
db1Ctr_Test = mean(db2Z_Test);

% Calculates the centroid to centroid coordinates
db1Proj_H0      = ((db2Z_H0 - db1Ctr_H0) * (db1Ctr_Test - db1Ctr_H0)' ./ ...
    norm(db1Ctr_Test - db1Ctr_H0).^2);
db1Proj_Test    = ((db2Z_Test - db1Ctr_H0) * (db1Ctr_Test - db1Ctr_H0)' ./ ...
    norm(db1Ctr_Test - db1Ctr_H0).^2);