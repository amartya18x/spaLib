% function [x_t t err] = tstage_lin(b, k, l, mxitr, tol, vtol)
function [x_t errNorm errSetDiff numItr] = tstage_lin(A, b, xTrue, k, l, mxitr, tol, vtol)
%Return the result of solving min \| b - Ax\|_2^2 , s.t support(x) <= k
%mxitr--maximum number of iterations, tol- tolerance for error=|b-Ax|
%vtol- tolerance for change in error
% global A;

[~, n]=size(A);
% eta=1;
% verbosity=0;

% err=+inf;
oerr=+inf;

x_all=A'*b;
[~, x_idx]=sort(abs(x_all),'descend');
curr_idx=x_idx(1:k);

x=A(:,curr_idx)\b;
%rp=randperm(n);
%curr_idx=rp(1:k);

x_t=sparse(n,1);

for t=1:mxitr
    
    res=A(:,curr_idx)*x-b;
    err=norm(res);
%     fprintf('%f\n',err);
    if err-oerr>tol
        break;
    end
    
    if err < tol || abs(oerr - err) < vtol
        break;
    end

%     [~, grad_idx]=GetMaxIdxLin(res, l);
    [~, grad_idx]=GetMaxIdxLin(A, res, l);
    new_idx=unique([grad_idx(:);curr_idx(:)]);
    z=A(:,new_idx)\b;
    [~, sortzidx]=sort(abs(z),'descend');
    curr_idx=new_idx(sortzidx(1:k));
    x=A(:,curr_idx)\b;

    oerr = err;
end
x_t(curr_idx)=x;

errNorm = norm(x_t - xTrue)/norm(xTrue);
errSetDiff = length(setdiff(find(xTrue),find(x_t)));
numItr = t;