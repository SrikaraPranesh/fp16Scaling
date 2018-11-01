% This script computes the change in the percentage of nonzero
% entries of a matrix because of underflow. Generates Table 4.5
% of the manuscript.

clear all; close all;
rng(1);

%%% Matlab file containing all the test matrices
load test_mat.mat



%%% prec_set=1 is (half,single,double) in GMRES-IR
%%% prec_set=2 is (half,double,quad  ) in GMRES-IR
prec_set = [1;2];

%%% theta -- Headroom to prevent overflow
theta = 0.1;


%%% flag=1 performs diagonal scaling in Algorithm 3.2
%%% flag=0 no diagonal scaling in Algorithm 3.2
flag = [0;1];

%%% dscale=1, uses Algorithm 2.4 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
%%% dscale=2, uses Algorithm 2.5 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
%%% dscale=3, uses Algorithm 2.6 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
dscale = [1;2;3];


%%%%% condtion number of diagonal scaling algorithms
for prec = 1:2
    for alg = 1:2
        for i = 1:length(test_mat)
            
            fprintf(' Matrix %d | Algorithm %d | Precisin %d \n',i,alg,prec);
            
            load(test_mat{i,1});  %%% Load the required matrix
            if (prec_set(prec,1) == 1)
                A = Problem.A;
                A = single(full(A));
                A = single(A);
                n = length(A);
                b = single(randn(n,1));
                u  =  eps('single');
                [u1,rmins,rmin,rmax,p] = ieee_params('h');
                rmax2 = single(rmax)*single(theta);
            else
                A = Problem.A;
                A = single(full(A));
                n = length(A);
                b = (randn(n,1));
                u = eps('double');
                [u1,rmins,rmin,rmax,p] = ieee_params('h');
                rmax2 = rmax*theta;
            end
            if (alg == 1)
                [B1,R1,C1] = scale_diag_2side(A);
                [B2,R2,C2] = scale_diag_2side_symm(A);
                
                B1 = double(fp16(B1));
                B2 = double(fp16(B2));
                
                cnnz{prec,alg}(i,1)=(abs(nnz(A) - nnz(B1))/nnz(A))*100;
                cnnz{prec,alg}(i,2)=(abs(nnz(A) - nnz(B2))/nnz(A))*100;
                
            elseif (alg == 2)
                
                [B1,R1,C1] = scale_diag_2side(A);
                mu1 = rmax/max(max(abs(B1)));
                [B2,R2,C2] = scale_diag_2side_symm(A);
                mu2 = rmax/max(max(abs(B2)));
                
                B1 = double(fp16(mu1*B1));
                B2 = double(fp16(mu2*B2));
                
                cnnz{prec,alg}(i,1)=(abs(nnz(A) - nnz(B1))/nnz(A))*100;
                cnnz{prec,alg}(i,2)=(abs(nnz(A) - nnz(B2))/nnz(A))*100;
                
            end

        end
        
    end
end


fid1=fopen('underflow_of_entries.txt','w');

fprintf(fid1,'algorithm 1 -- Un Symmetric diagonal scaling \n');
fprintf(fid1,'algorithm 2 -- Symmetric diagonal scaling \n');
fprintf(fid1,'\n'); fprintf(fid1,'\n');

for prec = 1:2
    fprintf(fid1,'Change in the percentage of non zeros in precision %d \n',alg,prec);
    for i=1:length(test_mat)
        t1 = cnnz{prec,1}(i,1); t2 = cnnz{prec,2}(i,1);
        t3 = cnnz{prec,1}(i,2); t4 = cnnz{prec,2}(i,2);
        fprintf(fid1,'%d & %.2f & %.2f & %.2f & %.2f\\\\ \n',i,...
            t1,t2,t3,t4);
    end
    fprintf(fid1,'\n'); fprintf(fid1,'\n');
end

fclose(fid1);

