function [A,mu] = scale_simple(A,rmax2,alg)
% This function performs the scaling as enumerated in
% Algorithm 2.1 and Algorithm 2.2
% A -- Input matrix
% rmax2 -- effective upper limit of fp16.
% alg=1 -- Conversion of the matrix using Algorithm 2.1
% alg=2 -- Conversion of the matrix using Algorithm 2.2

if (alg==1)    
    %%%% Conversion of the matrix using Algorithm 2.1
   j = find(abs(A) >= rmax2);
   A(j) = sign(A(j))*rmax2;
   mu = 1;
elseif (alg==2)
    %%%% Conversion of the matrix using Algorithm 2.2
    A1=zeros(length(A));
    amax = max(max(abs(A)));
    if (amax >= rmax2)
        mu = rmax2/amax;
    else
        mu=1;
    end
    A = mu*A;
end
end