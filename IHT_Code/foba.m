function [thetaFoBa errNorm errSetDiff numItr] = foba(X, y, thetaTrue, eps, expansion)
    nVars = length(thetaTrue);
    v = 0.5;
    supp = [];
    thetaFoBa = zeros(nVars,1);
    
    k = 0;
    numItr = 0;
    
    while length(supp) < expansion*nnz(thetaTrue)
        numItr = numItr+1;
        % Forward step
        currentLoss = getLossLasso(X,y,thetaFoBa);
        
        bestCoord = getBestCoord(X,y,thetaFoBa);
        if bestCoord < 0
            % No more coordinates to be added
            break;
        end
        
        % Add this newfound coordinate to the support set
        thetaFoBaTemp = thetaFoBa;
        supp = [supp, bestCoord];
        
        % Now apply fully corrective step
        thetaFoBa = solveRestrictedLS(X,y,supp);
        newLoss = getLossLasso(X,y,thetaFoBa);
        
        if currentLoss - newLoss <= eps
            % We made a mistake, remove this newly added coordinate and break out
            thetaFoBa = thetaFoBaTemp;
            break;
        else
            k = k+1;
            changeInLoss(k) = currentLoss - newLoss;
            currentLoss = newLoss;
        end	
        
        if length(supp) == expansion*nnz(thetaTrue)
            % No more edges to be added
            break;
        end
        
        % Backward steps
        while(true)
            worstCoord = getWorstCoord(X,y,thetaFoBa);
            if worstCoord < 0
                % No edges to be removed
                break;
            end
            
            thetaFoBaTemp = thetaFoBa;
            supp = supp(supp~=worstCoord);
            
            % Now apply fully corrective step
            thetaFoBa = solveRestrictedLS(X,y,supp);
            newLoss = getLossLasso(X,y,thetaFoBa);
            removeChange =  newLoss - currentLoss;
            %v=.1;
            if(removeChange > v*changeInLoss(k))
                % Removing this coordinate was a mistake. Add it back and break out to the outer loop
                thetaFoBa = thetaFoBaTemp;
                supp = [supp, worstCoord];
                break;
            else
                k = k-1;
                currentLoss = newLoss;
            end
        end
    end
    
    errNorm = norm(thetaFoBa - thetaTrue)/norm(thetaTrue);
    errSetDiff = length(setdiff(find(thetaTrue),find(thetaFoBa)));
end

function [loss] = getLossLasso(X,y,theta)
    loss = norm(X*theta - y)^2;
end

function [bestCoord] = getBestCoord(X,y,theta)
	bestCoord = -1;
    
    grad = getGradientLasso(X,y,theta);
    grad(theta ~= 0) = 0;
    
    if nnz(grad) == 0
        % We are out of coordinates to push in
        return;
    end
    
    [~, bestCoord] = max(abs(grad));
end

function [grad] = getGradientLasso(X,y,theta)
    grad = X'*(X*theta - y);
end

function [worstCoord] = getWorstCoord(X,y,theta)
    worstCoord = -1;
    res = y - X*theta;
    deltas = X(:,theta > 0);
    candidates = find(theta > 0);
    
    if isempty(deltas)
        return;
    end
    
    deltas = bsxfun(@times,theta(theta > 0)',deltas);
    newRes = bsxfun(@plus,deltas,res);
    newErr = sum(newRes.^2,1);
    [~, worstCoordLoc] = min(newErr);
    worstCoord = candidates(worstCoordLoc);
end

function [theta] = solveRestrictedLS(X,y,supp)
    thetaNewVals = X(:,supp)\y;
    theta = zeros(size(X,2),1);
    theta(supp) = thetaNewVals;
end
    