pBase = 20000;
sBase = 100;
oBase = 2;
eBase = 0.1;

numRuns = 3;
eta = 0.5;
mxitr = 100;
tol = 1e-3;
vtol = 1e-4;

% eps = [0.25 0.5 0.75 1 1.25 1.5 1.75 1.8];
eps = [1.75 1.8];
numDatasets = length(eps);
sizes = [1400];
numSizes = length(sizes);
sFactors = [0.75 1 1.25 1.5 1.75];
numSFactors = length(sFactors);

errNorm = zeros(numDatasets,numSizes,numSFactors,numRuns);
errDiff = zeros(numDatasets,numSizes,numSFactors,numRuns);
numItr = zeros(numDatasets,numSizes,numSFactors,numRuns);
recTime = zeros(numDatasets,numSizes,numSFactors,numRuns);

save('IHTErrNorm.mat','errNorm');
save('IHTErrDiff.mat','errDiff');
save('IHTNumItr.mat','numItr');
save('IHTRecTime.mat','recTime');

% Vary eps and generate data
for algoCtr = 1:numDatasets
    fprintf('HTP eps = %d\n',eps(algoCtr));
    filename = sprintf('Samples_p%d_s%d_e%g_o%g_Ceps%g.mat',pBase,sBase,eBase,oBase,eps(algoCtr));
    load(filename);
    for sizeCtr = 1:numSizes
        fprintf('\tsize %d\n',sizes(sizeCtr));
        for sCtr = 1:numSFactors
            fprintf('\t\tsFactor %d:',sFactors(sCtr));
            for run = 1:numRuns
                fprintf(' %d',run);
        
                tStart = tic;
                [temp, errNorm(algoCtr,sizeCtr,sCtr,run), errDiff(algoCtr,sizeCtr,sCtr,run) numItr(algoCtr,sizeCtr,sCtr,run)] = iiht_lin(samples.X(1:sizes(sizeCtr),:), samples.y(1:sizes(sizeCtr),run), samples.theta, sFactors(sCtr)*sBase, sFactors(sCtr)*sBase ,eta, mxitr, tol, vtol);
                recTime(algoCtr,sizeCtr,sCtr,run) = toc(tStart);

                save('IHTErrNorm.mat','errNorm');
                save('IHTErrDiff.mat','errDiff');
                save('IHTNumItr.mat','numItr');
                save('IHTRecTime.mat','recTime');
            end
            fprintf('\n');
        end
    end
end