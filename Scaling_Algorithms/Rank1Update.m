function [ C,b,alpha,beta,C1,R] = Rank1Update( A,b,theta,flag,dscale,rank1_type,mu)
% Performs diagonal scaling followed by rank-1 variant of
%scaling.
% A -- Input matrix.
% theta -- 0 < theta <= 1 headroom to prevent overflow.
% b - Right hand side vector
% flag - 1 if diagonal scaling is required.
% dscale - Type of diagonal scaling required
% rank1_type1 - Type of rank1 scaling
% mu - mu >=1 or 0 to modify the lower limit of scaled matrix.
% C - Input Matrix
% C1 - Column Scaling diagonal matrix
% R - Row scaling diagonal matrix 
% alpha, beta - Scaling constants


n = length(A);
%%% Equilibriate the matrix
if (flag == 1)
    if (dscale == 1)
        [A,R,C1] = scale_diag_2side(A);
    elseif (dscale == 2)
        [A,R,C1] = scale_diag_2side_symm(A);
    else
        [A,R,C1] = scale_diag_2side_symm_gm(A);
    end
else
    R = eye(n); C1 = eye(n);
end

[u1,rmins,rmin,rmax,p] = ieee_params('h');

if (rank1_type==1)
    rmax2 = rmax*theta;  % Allow some headroom.
    rmin2 = rmin*mu;  % Allow some headroom.
    AA = abs(A);
    xmax = max(AA(:));
    xmin = min(AA(:)) ;
    alpha = (rmax2 - rmin2)/(xmax - xmin);
    beta = (xmax*rmin2 - xmin*rmax2)/(xmax - xmin);
    C=alpha*A + beta*ones(n);  
else
    amax = max(A(:));
    amin = min(A(:));    
    alpha = (2*theta*rmax)/(amax-amin);
    beta = (theta*rmax*(-amin-amax))/(amax - amin);
    C = alpha*A + beta*ones(n);
end

b = R*b;

end

