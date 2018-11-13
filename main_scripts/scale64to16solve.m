function [x,gmresits,irits,Cnumber] = scale64to16solve(A,b,sf,maxit,theta,Scale,dscale)
%SCALE64TO16SOLVE  Solution of linear system of equation using scaling and
%GMRES based iterative refinement
%     x = scale64to16solve(A,b,sf,maxit,theta,Scale,dscale) Solves a system 
%     of linear equations using GMRES based iterative refinement. The input 
%     matrix is scaled to fit into the fp16 range, and the corresponding LU 
%     factors are used as preconditioners. 
%
%   Note: Requires Cleve Laboratory and Advanpix multiprecision toolbox

if (sf==1)
    A = single(A);
    b = single(b);
    tol = 1e-2;
    wp = 1;
    rp = 2;
    u = eps('single');
    [u1,rmins,rmin,rmax,p] = ieee_params('h');
    rmax2 = single(rmax)*single(theta);
else
    b = b;
    tol = 1e-4;
    wp = 2;
    rp = 4;
    u = eps('double');
    [u1,rmins,rmin,rmax,p] = ieee_params('h');
    rmax2 = rmax*theta;
end

if (Scale ==11 || Scale ==12)
    %%% Performs scaling using Algorithm 2.1 and Algorithm 2.2.
    if (Scale==11)
        [A1,mu] = SimpleScaling( A,rmax2,1 );
    elseif (Scale==12)
        [A1,mu] = SimpleScaling( A,rmax2,2);
    end
    Cnumber=cond(double(fp16(A1)),inf);
    if (Cnumber <= 1e20)
        [x,gmresits,irits] = gmresir3_1(double(A),double(b),0,wp,rp,maxit,tol,A1,mu);
    else
        x = zeros(length(b),1);
        gmresits = Inf;
        irits = Inf;
    end
    
elseif (Scale==2)
    %%% Performs scaling using just algorithm 2.3.
    Cnumber(1,1) = cond(A,inf);
    A_orig = A;
    b_orig = b;
    [ A,b,R,C,mu] = Diagonal_Scaling( A,b,dscale,rmax2 );
    Cnumber(1,2) = cond(A,inf);
    Cnumber(1,3) = cond(double(fp16(A)),inf);
    bInf = norm(b,'inf');
    b1 = b/bInf;
    if (Cnumber(1,3) <= (1e20) )
        [x,gmresits,irits]=gmresir3(double(A),double(b1),0,wp,rp,maxit,...
            tol,A_orig,b_orig,C,R,bInf,mu);
    else
        x = zeros(length(b),1);
        gmresits = Inf;
        irits = Inf;
    end
        
end

end

