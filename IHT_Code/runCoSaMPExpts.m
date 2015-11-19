pBase = 20000;
sBase = 100;
oBase = 2;
eBase = 0.1;

numDatasets = 8;
numRuns = 5;
mxitr = 100;
tol = 1e-3;
vtol = 1e-4;

errNorm = zeros(numDatasets,5);
errDiff = zeros(numDatasets,5);
numItr = zeros(numDatasets,5);
recTime = zeros(numDatasets,5);
thetaCoSaMP = zeros(pMax,numDatasets,5);

save('CoSaMPErrNorm.mat','errNorm');
save('CoSaMPErrDiff.mat','errDiff');
save('CoSaMPNumItr.mat','numItr');
save('CoSaMPRecTime.mat','recTime');
save('CoSaMPTheta.mat','thetaCoSaMP');

algoCtr = 0;

% Vary eps and generate data
for eps = [0.25 0.5 0.75 1 1.25 1.5 1.75 1.8]
    fprintf('CoSaMP eps =%d\n',eps);
    algoCtr = algoCtr + 1;
    filename = sprintf('Samples_p%d_s%d_e%g_o%g_Ceps%g.mat',pBase,sBase,eBase,oBase,eps);
    load(filename);
    for run = 1:numRuns
        fprintf('\tRun %d\n',run);
        sFactor = 1;
        if eps > 1.25
            sFactor = 2;
        end
        
        tStart = tic;
        [thetaCoSaMP(1:pBase,algoCtr,run), errNorm(algoCtr,run), errDiff(algoCtr,run) numItr(algoCtr,run)] = tstage_lin(samples.X, samples.y(:,run), samples.theta, sFactor*sBase, sFactor*2*sBase , mxitr, tol, vtol);
        recTime(algoCtr,run) = toc(tStart);
        
        save('CoSaMPErrNorm.mat','errNorm');
        save('CoSaMPErrDiff.mat','errDiff');
        save('CoSaMPNumItr.mat','numItr');
        save('CoSaMPRecTime.mat','recTime');
    end
end
save('CoSaMPTheta.mat','thetaCoSaMP');