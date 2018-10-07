function [ A,b,R,C ] = Diagonal_Scaling( A,b,dscale,rmax2 )
% This function performs the scaling as enumerate in
% algorithm 2.3.
% A - The input matrix
%%% Performs diagonal scaling using either
%%% Algorithm 2.4 or 2.5.
%%% dscale = 1 -- Algorithm 2.4
%%% dscale = 2 -- Algorithm 2.5
%%% dscale = 3 -- Algorithm 2.6

if (dscale==1)
    [A,R,C] = scale_diag_2side(A);
elseif (dscale==2)
    [A,R,C] = scale_diag_2side_symm(A);
elseif (dscale==3)
    [A,R,C] = scale_diag_2side_symm_gm(A);
end
beta=max(max(A));
mu=rmax2/beta;
A=mu*A;
b=mu*R*b;

end

