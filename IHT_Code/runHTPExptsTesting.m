pBase = 20000;
sBase = 100;
oBase = 2;
eBase = 0.1;

numDatasets = 7;
numRuns = 5;
eta = 0.5;
mxitr = 100;
tol = 1e-3;
vtol = 1e-4;

errNorm = zeros(numDatasets,5);
errDiff = zeros(numDatasets,5);
numItr = zeros(numDatasets,5);
recTime = zeros(numDatasets,5);
bestSFactors = zeros(numDatasets,5);
thetaIHT = zeros(pBase,numDatasets,5);

save('IHTErrNorm.mat','errNorm');
save('IHTErrDiff.mat','errDiff');
save('IHTNumItr.mat','numItr');
save('IHTRecTime.mat','recTime');
save('IHTSFactor.mat','bestSFactors');
save('IHTTheta.mat','thetaIHT');

algoCtr = 0;

% Vary e and generate data
for eps = [0.25 0.5 0.75 1 1.25 1.5 1.75]
    fprintf('HTP eps=%d\n',eps);
    algoCtr = algoCtr + 1;
    filename = sprintf('Samples_p%d_s%d_e%g_o%g_Ceps%g.mat',pBase,sBase,eBase,oBase,eps);
    load(filename);
    for run = 1:numRuns
        fprintf('\tRun %d\n',run);
        errDiffBest = inf;
        sFactorBest = inf;
        thetaIHTBest = [];
        erroNormBest = inf;
        numItrBest = inf;
        recTimeBest = inf;
        
        for sFactor = [1 2 4 8]
            fprintf('\t\tsFactor %d\n',sFactor);
            tStart = tic;
            [thetaIHTTemp, errNormTemp, errDiffTemp, numItrTemp] = iiht_lin(samples.X, samples.y(:,run), samples.theta, sFactor*sBase, sFactor*sBase ,eta, mxitr, tol, vtol);
            recTimeTemp = toc(tStart);
            
            if errDiffTemp < errDiffBest
                errDiffBest = errDiffTemp;
                thetaIHTBest = thetaIHTTemp;
                erroNormBest = errNormTemp;
                numItrBest = numItrTemp;
                recTimeBest = recTimeTemp;
                sFactorBest = sFactor;
            end
        end
        
        thetaIHT(:,algoCtr,run) = thetaIHTBest;
        errNorm(algoCtr,run) = erroNormBest;
        errDiff(algoCtr,run) = errDiffBest;
        numItr(algoCtr,run) = numItrBest;
        recTime(algoCtr,run) = recTimeBest;
        bestSFactors(algoCtr,run) = sFactorBest;
        
        save('IHTErrNorm.mat','errNorm');
        save('IHTErrDiff.mat','errDiff');
        save('IHTNumItr.mat','numItr');
        save('IHTRecTime.mat','recTime');
        save('IHTSFactor.mat','bestSFactors');
    end
end
save('IHTTheta.mat','thetaIHT');