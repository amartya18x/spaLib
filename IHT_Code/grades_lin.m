% function [x t err] = grades_lin(b, k, eta,mxitr, tol, vtol)
function [x errNorm errSetDiff numItr] = grades_lin(A, b, xTrue, k, eta, mxitr, tol, vtol)
%Return the result of solving min \| b - Ax\|_2^2 , s.t support(x) <= k
%mxitr--maximum number of iterations, tol- tolerance for error=|b-Ax|
%vtol- tolerance for change in error
% global A;
[~, n] = size(A);
x=zeros(n,1);
verbosity=0;

% err=+inf;
oerr=+inf;
if verbosity == 1
    fprintf('\nIteration:   ');
end

for t=1:mxitr
    
    x=x-eta*A'*(A*x-b);
    [~, sortidx]=sort(abs(x),'descend');
    x(sortidx(k+1:n))=0;
    err=norm(A*x - b);
    %incorr_idx=length(setdiff(sortidx(1:k), find(x_opt)));
    if verbosity == 1
        fprintf('\b\b\b\b%4d',t);
    elseif verbosity >1
        fprintf('Iteration %d: Error = %f\n', t, err);
    end
    if err-oerr>tol
        err=norm(b);
        break;
    end
    if err < tol || abs(oerr - err) < vtol
        break;
    end

    oerr = err;
end
% fprintf('\n');

errNorm = norm(x - xTrue)/norm(xTrue);
errSetDiff = length(setdiff(find(xTrue),find(x)));
numItr = t;