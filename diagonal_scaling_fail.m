% This script computes the condition number of the scaled matrix for
% varying values of \mu. Results will be printed into 
% 'IntervalVary.txt'.


clear all
close all


addpath('matrices_for_testing')
addpath('GMRES_IR')
addpath('DiagScaling')
addpath('Scaling_Algorithms')
addpath('Misc_Scripts')


load test_mat_2.mat

for i = 1:3
    load(test_mat_2{i,1});
    A = Problem.A;
    A = A((1:300),(1:300));
    A = full(A);
    dscale = 1;
    [u1,rmins,rmin,rmax,p] = ieee_params('h');
    rmax2 = rmax*0.1;  % Allow some headroom.
    
    %%%% diagonal scaling
    if (dscale==1)
        [B,R,C1] = scale_diag_2side(A);
    elseif (dscale==2)
        [B,R,C1] = scale_diag_2side_symm(A);
    end
    
    beta = max(max(B));
    mu = rmax2/beta;
    B1 = B;
    B = mu*B;
    Bh = double(fp16(B));
    
    val = [10;10^5];
    
    %%%% rank1 scaling
    for j=1:2
        rmin2 = rmin*val(j,1);  % Allow some headroom.
        AA = abs(B1);
        xmax = max(AA(:));
        xmin = min(AA(:)) ;
        alpha = (rmax2 - rmin2)/(xmax - xmin);
        beta = (xmax*rmin2 - xmin*rmax2)/(xmax - xmin);
        C = alpha*B1 + beta*ones(length(B1));
        Ch = double(fp16(C));
        
        rval{i,j} = [cond(inv(Bh)*B,inf) cond(inv(Ch)*C,inf);
            cond(B,inf) cond(Bh,inf);
            cond(C,inf) cond(Ch,inf)];
    end
    
end


fid1 = fopen('IntervalVary.txt','w');

%%%%% Print condition number of scaled matrix into a file
fprintf(fid1,'Condition number of the scaled matrix for varying values of \mu')
fprintf(fid1,'\\mu=0 & \\mu=10 \\times \\xmin & \\mu=10^5 \\times \\xmin \n')
for i = 1:3
    fprintf(fid1,'rval{i,1}(2,1) & rval{i,1}(3,1) & rval{i,2}(3,1) \\\\ \n');
end
fprintf(fid1,'\n');fprintf(fid1,'\n');

%%%%% Print condition number of scaled fp16 matrix into a file
fprintf(fid1,'Condition number of the scaled fp16 matrix for varying \mu')
fprintf(fid1,'\\mu=0 & \\mu=10 \\times \\xmin & \\mu=10^5 \\times \\xmin \n')
for i = 1:3
    fprintf(fid1,'rval{i,1}(2,2) & rval{i,1}(3,2) & rval{i,2}(3,2) \\\\ \n');
end
fprintf(fid1,'\n');fprintf(fid1,'\n');

%%%%% Print condition number of scaled fp16 matrix into a file
fprintf(fid1,'Condition number of the left preconditioned matrix for varying \mu')
fprintf(fid1,'\\mu=0 & \\mu=10 \\times \\xmin & \\mu=10^5 \\times \\xmin \n')
for i = 1:3
    fprintf(fid1,'rval{i,1}(1,1) & rval{i,1}(1,2) & rval{i,2}(1,2) \\\\ \n');
end





