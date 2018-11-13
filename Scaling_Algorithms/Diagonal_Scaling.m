function [ A,b,R,C,mu ] = Diagonal_Scaling( A,b,dscale,rmax2 )
%DIAGONAL_SCALING Performs diagonal scaling and multiplies by 
%(theta_max X xmax)

global mu_flag;

if (dscale==1)
    [A,R,C] = scale_diag_2side(A);
elseif (dscale==2)
    [A,R,C] = scale_diag_2side_symm(A);
elseif (dscale==3)
    [A,R,C] = scale_diag_2side_symm_gm(A);
end
beta = max(max(A));
mu = rmax2/beta;

if (mu_flag == 1)
   mu = 1; 
end
A = mu*A;
b = mu*R*b;

end

