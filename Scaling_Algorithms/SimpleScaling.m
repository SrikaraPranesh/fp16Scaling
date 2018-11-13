function [ A1,mu ] = SimpleScaling( A,rmax2,alg )
% SIMPLESCALING Scaling by mapping to maximum number, or divide by
% maximum number.

if (alg==1)    
    %%%% Conversion of the matrix using Algorithm 2.1
    A1 = zeros(length(A));
    for i1 = 1:length(A)
        for j = 1:length(A)
            if (abs(A(i1,j)) >= rmax2)
                A1(i1,j) = rmax2;
            else
                A1(i1,j) = A(i1,j);
            end
        end
    end
    mu = 1;
elseif (alg==2)
    %%%% Conversion of the matrix using Algorithm 2.2
    A1 = zeros(length(A));
    amax = max(max(abs(A)));
    if (amax >= rmax2)
        mu = rmax2/amax;
    else
        mu = 1;
    end
    A1 = mu*A;
end




end

