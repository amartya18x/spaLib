% function [val idx]=GetMaxIdxLin(query,k)
function [val idx]=GetMaxIdxLin(A,query,k)
% global A;
[val idx]=sort(abs(A'*query),'descend');
%toc;
val=val(1:k);
idx=idx(1:k);

