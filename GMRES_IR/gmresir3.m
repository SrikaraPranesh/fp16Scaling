function [x1,reqits,iter] = gmresir3(A,b,precf,precw,precr,iter_max,gtol,A1,b1,C,R,bInf,mu)
%GMRESIR3  GMRES-based iterative refinement in three precisions.
%     x = gmresir3(A,b,precf,precw,precr,iter_max,gtol) solves Ax = b using gmres-based
%     iterative refinement with at most iter_max ref. steps and GMRES convergence
%     tolerance gtol, with
%     LU factors computed in precision precf:
%       * half if precf = 0,
%       * single if precf = 1,
%       * double if precf = 2,
%     working precision precw:
%       * half if precw = 0,
%       * single if precw = 1,
%       * double if precw = 2,
%     and residuals computed at precision precr:
%       * single if precr = 1,
%       * double if precr = 2,
%       * quad if precr = 4
%
%   Note: Requires Cleve Laboratory and Advanpix multiprecision toolbox

if precf ~=0 && precf ~=1 && precf ~= 2, error('precf should be 0, 1 or 2'), end
if precw ~=0 && precw ~=1 && precw ~= 2, error('precw should be 0, 1 or 2'), end
if precr ~=1 && precr ~= 2 && precr ~= 4, error('precr should be 1, 2, or 4'), end

n = length(A);

if precf == 1
    fprintf('**** Factorization precision is single.\n')
    ufs = 'single';
elseif precf == 2
    fprintf('**** Factorization precision is double.\n')
    ufs = 'double';
else
    fprintf('**** Factorization precision is half.\n')
    ufs = 'half';
end

if precw == 0
    fprintf('**** Working precision is half.\n')
    uws = 'half';
    A = fp16(A);
    b = fp16(b);
    u = eps(fp16(1));
elseif precw == 2
    fprintf('**** Working precision is double.\n')
    uws = 'double';
    A = double(A);
    b = double(b);
    u = eps('double');
else
    fprintf('**** Working precision is single.\n')
    uws = 'single';
    A = single(A);
    b = single(b);
    u = eps('single');
end

if precr == 1
    fprintf('**** Residual precision is single.\n')
    urs = 'single';
elseif precr == 2
    fprintf('**** Residual precision is double.\n')
    urs = 'double';
else
    fprintf('**** Residual precision is quad.\n')
    urs = 'quad';
    mp.Digits(34);
end

xact = double(mp(double(A1),34)\mp(double(b1),34));

%Compute LU factorization
if precf == 1
    [L,U,P] = lu(single(A));
    LL = single(double(P')*double(L));
    x =  U\(L\(P*single(b)) );
elseif precf == 2
    [L,U,P] = lu(double(A));
    LL = double(double(P')*double(L));
    x =  U\(L\(P*double(b)) );
else
    B=double(fp16(A));
    [L,U,P] = lu(B);
    LL = fp16(double(P')*double(L));
    U=fp16(U);
    x =  U\(L\(P*fp16(b)) );
    
    LL = (1/mu)*diag(1./diag(R))*double(LL);
    U = double(U)*diag(1./diag(C));
end

%Compute condition number of A, of preconditioned system At, cond(A), and
%cond(A,x) for the exact solution
At = double(mp(double(U),34)\(mp(double(L),34)\( mp(double(P),34)*mp(double(A),34))));
kinfA = cond(mp(double(A),34),'inf');
kinfAt = cond(mp(double(At),34),'inf');

%Note: when kinf(A) is large, the initial solution x can have 'Inf's in it
%If so, default to using 0 as initial solution
if sum(isinf(single(x)))>0
    x =  zeros(size(b,1),1);
    fprintf('**** Warning: x0 contains Inf. Using 0 vector as initial solution.\n')
end

%Store initial solution in working precision
if precw == 0
    x = fp16(x);
    x1=bInf*C*x;
elseif precw == 2
    x = double(x);
    x1=bInf*C*x;
else
    x = single(x);
    x1=bInf*C*x;
end

cged = false;
iter = 0; dx = 0; rd = 0; sflag=0;

%Array to store total number of gmres iterations in each ref step
gmresits = [];
reqits=0;
%Array to store final relative (preconditioned) residual norm in gmres
gmreserr = [];

while ~cged
    
    %Compute size of errors, quantities in bounds
    ferr(iter+1) = double(norm(mp(double(x1),34)-mp(xact,34),'inf')/norm(mp(xact,34),'inf'));
    res = double(b1) - double(A1)*double(x1);
    nbe(iter+1) = double(norm(mp(res,34),'inf')/(norm(mp(double(A1),34),'inf')*norm(mp(double(x1),34),'inf')+ norm(mp(double(b1),34),'inf')));
    
    iter = iter + 1;
    if iter > iter_max, break, end
    
    %Check convergence
%     if max([ferr(iter) nbe(iter)]) <= (length(b1)*u), break, end
    if (([ nbe(iter) ]) <= (length(b1)*u)), break, end
    
    %Compute residual vector
    if precr == 1
        rd = single(b1) - single(A1)*single(x1);
    elseif precr == 2
        rd = double(b1) - double(A1)*double(x1);
    else
        rd = mp(double(b1),34) - mp(double(A1),34)*mp(double(x1),34);
    end
    
    %Scale residual vector
    norm_rd = norm(rd,inf);
    rd1 = rd/norm_rd;
    A11 = A1;
    %Call GMRES to solve for correction term
    if precw == 0
        [d, err, its, ~] = gmres_hs( A11, fp16(zeros(n,1)), fp16(rd1), LL, U, n, 1, gtol);
    elseif precw == 2
        [d, err, its, ~] = gmres_dq( A11, zeros(n,1), double(rd1), LL, U, n, 1, gtol);
    else
        [d, err, its, ~] = gmres_sd( A11, single(zeros(n,1)), single(rd1), LL, U, n, 1, gtol);
    end

    
    %Record number of iterations gmres took
    gmresits = [gmresits,its];
    reqits=reqits+its;
    %Record final relative (preconditioned) residual norm in GMRES
    gmreserr = [gmreserr,err(end)];
    
    %Record relative (preconditioned) residual norm in each iteration of
    %GMRES (so we can look at convergence trajectories if need be)
    gmreserrvec{iter} = err;
    
    xold = x1;
    
    %Update solution
    if precw == 0
        x1 = x1 + fp16(norm_rd)*fp16(d);
    elseif precw == 2
        x1 = x1 + norm_rd*double(d);
    else
        x1 = x1 + single(norm_rd)*single(d);
    end
    dx = norm(x1-xold,'inf')/norm(x1,'inf');
    
    %Check if dx contains infs, nans, or is 0
    if dx == Inf || isnan(double(dx))
        plt = 0;
        break;
    end
    
end

% if ((iter >= iter_max) && (max([nbe(iter) ferr(iter)]) > length(b)*u))
if ((iter >= iter_max) && (nbe(iter) > length(b)*u))
    reqits=Inf;
end

end
