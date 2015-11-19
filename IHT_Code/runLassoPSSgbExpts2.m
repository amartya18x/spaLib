pBase = 20000;
sBase = 100;
oBase = 2;
eBase = 0.1;

numDatasets = 18;
numRuns = 5;
eta = 0.5;
mxitr = 100;
tol = 1e-3;
vtol = 1e-4;

errNorm = zeros(numDatasets,5);
errDiff = zeros(numDatasets,5);
numItr = zeros(numDatasets,5);
recTime = zeros(numDatasets,5);
thetaLassoPSSgb = zeros(pBase,numDatasets,5);

save('LassoPSSgbErrNorm_2.mat','errNorm');
save('LassoPSSgbErrDiff_2.mat','errDiff');
save('LassoPSSgbNumItr_2.mat','numItr');
save('LassoPSSgbRecTime_2.mat','recTime');
save('LassoPSSgbTheta_2.mat','thetaLassoPSSgb');

algoCtr = 0;

% Vary eps and generate data
for eps = [0.25 0.5 0.75 1 1.25 1.5 1.75 1.8]
    fprintf('LassoPSSgb eps =%d\n',eps);
    algoCtr = algoCtr + 1;
    filename = sprintf('Samples_p%d_s%d_e%g_o%g_Ceps%g.mat',pBase,sBase,eBase,oBase,eps);
    load(filename);
    for run = 2:numRuns
        fprintf('\tRun %d\n',run);
        sFactor = 1;
        if numRuns > 5
            sFactor = 2;
        end
        
        [thetaLassoPSSgb(1:pBase,algoCtr,run), errNorm(algoCtr,run), errDiff(algoCtr,run) recTime(algoCtr,run) numItr(algoCtr,run)] = runLasso_PSSgb(samples.X, samples.y(:,run), samples.theta, 1/sqrt(size(samples.X,1)), sFactor);
        
        save('LassoPSSgbErrNorm_2.mat','errNorm');
        save('LassoPSSgbErrDiff_2.mat','errDiff');
        save('LassoPSSgbNumItr_2.mat','numItr');
        save('LassoPSSgbRecTime_2.mat','recTime');
    end
end
save('LassoPSSgbTheta_2.mat','thetaLassoPSSgb');