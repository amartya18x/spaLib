function [bestThetaLasso bestErrNorm bestErrSetDiff recoveryTimeForBestErrSetDiff bestNumItr] = runLasso_PSSgb(X, y, thetaTrue, lambda, expansion)
    bestErrSetDiff = inf;
    bestErrNorm = inf;
    recoveryTimeForBestErrSetDiff = inf;
    bestThetaLasso = [];
    
    [~, nVars] = size(X);
    funObj = @(theta)SquaredError(theta,X,y);
    
    for i = [1 2 5 10 20 40]
        tStart = tic;
        
        % Penalize the absolute value of each element by the same amount
        lambdaVec = i*lambda*ones(nVars,1);
        
        % Perform ridge regression to obtain a warm start value
        R = chol(X'*X + diag(lambdaVec));
        thetaRR = R\(R'\(X'*y));
        theta_init = thetaRR;
        
        % Get the LASSO estimate
        [thetaLasso numItr] = L1General2_PSSgb(funObj,theta_init,lambdaVec);
        
        % Get the top most coordinates of the LASSO estimate
        [~, mx_idx]=sort(abs(thetaLasso), 'descend');
        
        % Calculate recovery properties using this
        errSetDiff = length(setdiff(find(thetaTrue),mx_idx(1:expansion*nnz(thetaTrue))));
        thetaProj = thetaLasso;
        thetaProj(mx_idx(expansion*nnz(thetaTrue)+1:end)) = 0;
        errNorm = norm(thetaProj - thetaTrue)/norm(thetaTrue);
        
        recoveryTime = toc(tStart);
        
        % Get the best time possible for the best recovery
        if errSetDiff <= bestErrSetDiff
            if recoveryTime < recoveryTimeForBestErrSetDiff
                bestErrSetDiff = errSetDiff;
                bestErrNorm = errNorm;
                recoveryTimeForBestErrSetDiff = recoveryTime;
                bestNumItr = numItr;
                bestThetaLasso = thetaLasso;
            end
        end
    end
end

function [error grad] = SquaredError(theta,X,y)
    error = norm(X*theta - y)^2;
    grad = X'*(X*theta - y);
end