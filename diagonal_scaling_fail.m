clear all
close all


addpath('matrices_for_testing')
addpath('GMRES_IR')
addpath('DiagScaling')
addpath('Scaling_Algorithms')
addpath('Misc_Scripts')
addpath('../../Matrix_Collection/SuiteSparse')

% A=[1+(3*eps) 1+eps; 1-eps 1+(2*eps)];
%
load test_mat_2.mat

for i=1:3
    load(test_mat{i,1});
    A=Problem.A;
    A=full(A);
    dscale=1;
    [u1,rmins,rmin,rmax,p] = ieee_params('h');
    rmax2 = rmax*0.1;  % Allow some headroom.
    
    %%%% diagonal scaling
    if (dscale==1)
        [B,R,C1] = scale_diag_2side(A);
    elseif (dscale==2)
        [B,R,C1] = scale_diag_2side_symm(A);
    end
    
    beta=max(max(B));
    mu=rmax2/beta;
    B1=B;
    B=mu*B;
    Bh=double(fp16(B));
    
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
        
        rval{j,1} = [cond(inv(Bh)*B,inf) cond(inv(Ch)*C,inf);
            cond(B,inf) cond(Bh,inf);
            cond(C,inf) cond(Ch,inf)];
    end
    
end


fid1=fopen('IntervalVary.txt','w');


fprintf(fid1,'\mu= 10 \times \xmin\n');
fprintf(fid1,'\mu= 10 \times \xmin\n');


