function [ A1 ] = SimpleScaling( A,rmax2,alg )
% This function performs the scaling as enumerated in
% Algorithm 2.1 and Algorithm 2.2
% A -- Input matrix
% rmax2 -- effective upper limit of fp16.
% alg=1 -- Conversion of the matrix using Algorithm 2.1
% alg=2 -- Conversion of the matrix using Algorithm 2.2

if (alg==1)    
    %%%% Conversion of the matrix using Algorithm 2.1
    A1=zeros(length(A));
    for i1=1:length(A)
        for j=1:length(A)
            if (abs(A(i1,j)) >= rmax2)
                A1(i1,j)=rmax2;
            else
                A1(i1,j)=A(i1,j);
            end
        end
    end
elseif (alg==2)
    %%%% Conversion of the matrix using Algorithm 2.2
    A1=zeros(length(A));
    amax=max(max(abs(A)));
    if (amax >= rmax2)
        mu=rmax2/amax;
    else
        mu=1;
    end
    A1=mu*A;
end




end

