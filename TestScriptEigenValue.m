% Small test for eigen value expected value.
inNIter = 1000;
inNObs  = 150;
inNVar  = 50;

[db1Av_Eig, db1SD_Eig] = deal(nan(1, inNIter));
for iItr = 1:inNIter
    db2RandNorm = randn(inNObs, inNVar);
    db1Eig      = eig(db2RandNorm' * db2RandNorm);
    db1Av_Eig(iItr) = mean(db1Eig);
    db1SD_Eig(iItr) = std(db1Eig);
end
fprintf('Mean Mean:\t%.3f\rMean S.D.:\t%.3f\r', mean(db1Av_Eig), mean(db1SD_Eig))