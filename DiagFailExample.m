% An example where just the diagonal scaling
% will make the fp16 matrix singular.
% Will display the eigenvalues of scaled matrix,
% and one can observe that it is singular

clear all
close all


addpath('matrices_for_testing')
addpath('GMRES_IR')
addpath('DiagScaling')
addpath('Scaling_Algorithms')
addpath('Misc_Scripts')


A = [1+(3*eps) 1+eps; 1-eps 1+(2*eps)];
[u1,rmins,rmin,rmax,p] = ieee_params('h');
rmax2 = rmax*0.1;
%%%% diagonal scaling using equilibriation
[B1,R,C1] = scale_diag_2side(A);


beta = max(max(B1));
mu = rmax2/beta;
B1 = mu*B1;
Bh1 = double(fp16(B1));
eig(Bh1)

%%%% Diagonal scaling using symmetric scaling
[B2,R,C1] = scale_diag_2side_symm(A);



beta = max(max(B2));
mu = rmax2/beta;
B2 = mu*B2;
Bh2 = double(fp16(B2));
eig(Bh2)


