function db1MahaD = MahalanobisD(db2RefPop, db2Obs)
%DB1MAHAD = MahalanobisD(DB2REFPOP, DB2OBS)
%Returns DB1MAHAD: the Mahalanobis distances of a set of observation
%DB2OBS relative to a set of reference observation DB2REFPOP.

db2CCov     = nancov(db2RefPop); %Covariance matrix of the cluster
db1CMu      = nanmean(db2RefPop, 1); %Centroid of the cluster

db2ObsCnt   = db2Obs - db1CMu; %Centered coordintates
db1MahaD    = sqrt(nansum(db2ObsCnt/db2CCov.*db2ObsCnt, 2)); %Squared Mahalonobis distances