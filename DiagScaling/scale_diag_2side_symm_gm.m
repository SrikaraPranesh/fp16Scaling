function [A,D1,D2,its] = scale_diag_2side_symm_gm(A,tol,prnt)
%SCALE_DIAG_2SIDE_SYMM_GM Symmetry-preserving two-sided GM diagonal scaling.
%   [B,d1,d2,its] = scale_diag_2side_symm_gm(A,TOL) computes B = D1*A_D2,
%   where the diagonal matrices D1 = inv(diag(d1)) and D2 = inv(diag(d2))
%   are built from approximations to the geoemtic means of rows and
%   columns.  The scaling is symmetry-preserving.
%   A must not have a zero row or column.
%   TOL is a convergence tolerance: default TOL = 1e-4.
%   ITS is the number of iterations for convergence.

if nargin < 2 || isempty(tol), tol = 1e-4; end
if nargin < 3, prnt = 0; end

n = length(A);

d1prod = ones(n,1);
d2prod = ones(1,n);

A_old = zeros(n);

if prnt
   fprintf('%2.0f:  dA = 0,         normA = %9.2e, condA = %9.2e\n',...
           0, norm(A,1), cond(A,1))
end   
cc=1;
for k = 1:1

    A_old = A;
    absA_nozeros = abs(A) + realmax*(A == 0);
    row_min = min(absA_nozeros,[],2); row_max = max(abs(A),[],2);
    d1 = sqrt(row_min .* row_max);
    %A = diag(1./d1)*A;  % 
    col_min = min(absA_nozeros,[],1); col_max = max(abs(A),[],1);
    d2 = sqrt(col_min .* col_max);
    A = diag(1./d1)*A*diag(1./d2);  % Implicit expansion: same as A/diag(d2).
    d1prod = d1prod./d1;
    d2prod = d2prod./d2;
    if prnt
       fprintf('%2.0f:  dA = %9.2e, normA = %9.2e, condA = %9.2e, ConvergCrit= %9.2e\n',...
               k, norm(A-A_old,1), norm(A,1), cond(A,1),cc)
    end   
    cc = norm([row_max' col_max] - 1, inf);
    if cc < tol, its = k; break, end
    
end
D1=diag(d1prod); D2=diag(d2prod);
its = k;