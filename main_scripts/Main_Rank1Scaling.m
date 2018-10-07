% This code performs scaling of the scaling of a general
% matrix in double or single precision to half precision.
% The algorithm numbers used are consisten with the ones
% used in the paper "Squeezing a Matrix Into Half Precision,
% with an Application to Solving Linear Systems -- N.J Higham,
% S. Pranesh, and M. Zounon"

% This sybroutine, requires the following
% 1.  GMRES based iterative refinement codes, which can be
% found in https://github.com/eccarson/ir3
%
% 2. Cleve's Lab for fp16, and advanpix toolbox for
% high precision computation.
%
clear all;close all;


rank1_type = 1;
addpath('matrices_for_testing')
addpath('GMRES_IR')
addpath('DiagScaling')
addpath('Scaling_Algorithms')
addpath('Misc_Scripts')
rng(1);

%%% Matlab file containing all the test matrices
load test_mat.mat



%%% prec_set = 1 is (half,single,double) in GMRES-IR
%%% prec_set = 2 is (half,double,quad  ) in GMRES-IR
prec_set = [1;2];


%%% theta -- Headroom to prevent overflow
%%% mu -- Room to prevent underflow
theta = 0.1;
mu = 1;

%%% scale = 11 Performs scaling using Algorithm 2.1
%%% scale = 12 Performs scaling using Algorithm 2.2.
%%% scale = 2  Performs scaling using just algorithm 2.3.
%%% scale = 3 Performs scaling using algorithm 3.2
Scale = 3 ;


%%% flag = 1 performs diagonal scaling in Algorithm 3.2
%%% flag = 0 no diagonal scaling in Algorithm 3.2
flag = [0;1]  ;


%%% Total number iterative refinement iterations
maxit = 10;

%%% dscale = 1, uses Algorithm 2.4 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
%%% dscale = 2, uses Algorithm 2.5 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
%%% dscale = 3, uses Algorithm 2.6 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
dscale = 2;

for prec = 1:2
    
    for alg = 1:2
        
        for i = 1:length(test_mat)
            load(test_mat{i,1});  %%% Load the required matrix
            
            if (prec_set(prec,1)==1)
                A = Problem.A;
                A = single(full(A));
                n = length(A);
                b = single(randn(n,1));
                lm(i,1) = n;
            else
                A = Problem.A;
                A = (full(A));
                n = length(A);
                b = (randn(n,1));
                lm(i,1) = n;
            end
            
            fprintf('I am in matrix number %d \n',i);
            [x,gmresits{prec,alg}(i,:),irits{prec,alg}(i,:),Cnumber{prec,alg}(i,:)...
                tau{prec,alg}(i,:),absx{prec,alg}(i,:)] = ...
                Scale64To16_Solve(A,b,prec_set(prec,1)...
                ,maxit,theta,Scale,flag(alg,1),dscale,mu);
            res = double(b) - double(A)*double(x);
            nbe{prec,alg}(i,1) = double(norm(mp(res,34),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34),'inf')+ norm(mp(double(b),34),'inf')));
        end    

    end
    
end


%%%%%%% PRINT THE LATEX TABLE INTO A FILE %%%%%%%

% creating a text file to print the GMRES iteration table
if (rank1_type==1)
    fid1 = fopen('Rank1ScalingResult_type11.txt','w');
else
    fid1 = fopen('Rank1ScalingResult_type2.txt','w');
end
fprintf(fid1,'Number of GMRES iterations for (half,single,double) precisions combination, without diagonal scaling \n');
fprintf(fid1,'Index & $By=d$ & $Bz=e$ \\\\ \n')
for i = 1:length(test_mat)
    t1 = gmresits{1,1}(i,1); t2 = gmresits{1,1}(i,2);
    fprintf(fid1,'%d & %d & %d \\\\ \n',i,...
        t1,t2);
end
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'Number of GMRES iterations for (half,double,quad) precisions combination, without diagonal scaling \n');
fprintf(fid1,'Index & $By=d$ & $Bz=e$ \\\\ \n')
for i = 1:length(test_mat)
    t1 = gmresits{2,1}(i,1); t2 = gmresits{2,1}(i,2);
    fprintf(fid1,'%d & %d & %d \\\\ \n',i,...
        t1,t2);
end
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'Number of GMRES iterations for (half,single,double) precisions combination, with diagonal scaling \n');
fprintf(fid1,'Index & $By=d$ & $Bz=e$ \\\\ \n')
for i = 1:length(test_mat)
    t1 = gmresits{1,2}(i,1); t2 = gmresits{1,2}(i,2);
    fprintf(fid1,'%d & %d & %d \\\\ \n',i,...
        t1,t2);
end
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'Number of GMRES iterations for (half,double,quad) precisions combination, with diagonal scaling \n');
fprintf(fid1,'Index & $By=d$ & $Bz=e$ \\\\ \n')
for i = 1:length(test_mat)
    t1 = gmresits{2,2}(i,1); t2 = gmresits{2,2}(i,2);
    fprintf(fid1,'%d & %d & %d \\\\ \n',i,...
        t1,t2);
end
fprintf(fid1,'\n');
fprintf(fid1,'\n');

%%%%% printing backward errors
fprintf(fid1,'\n');
fprintf(fid1,'\n');
fprintf(fid1,'backward errors for Single precision \n');
fprintf(fid1,'Index & No diag scaling & Diag scaling \\\\ \n')
for i = 1:length(test_mat)
    t1 = nbe{1,1}(i,1); t2 = nbe{1,2}(i,1);
    fprintf(fid1,'%d & %6.2e & %6.2e \\\\ \n',i,...
        t1,t2);
end

fprintf(fid1,'\n');
fprintf(fid1,'\n');
fprintf(fid1,'backward errors double precision\n');
fprintf(fid1,'Index & No diag scaling & Diag scaling \\\\ \n')
for i = 1:length(test_mat)
    t1 = nbe{2,1}(i,1); t2 = nbe{2,2}(i,1);
    fprintf(fid1,'%d & %6.2e & %6.2e\\\\ \n',i,...
        t1,t2);
end


%%%%% terms in error analysis
fprintf(fid1,'\n');
fprintf(fid1,'\n');
fprintf(fid1,'Terms from error analysis in single precision without diagonal scaling\n');
fprintf(fid1,'Index & \tau & |h| & \frac{xprime}{xtilde} \\\\ \n')
for i = 1:length(test_mat)
    t1 = tau{1,1}(i,1); t2 = tau{1,1}(i,2); t3 = absx{1,1}(i,1);
    fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e\\\\ \n',i,...
        t1,t2,t3);
end

fprintf(fid1,'\n');
fprintf(fid1,'\n');
fprintf(fid1,'Terms from error analysis in double precision without diagonal scaling\n');
fprintf(fid1,'Index & \tau & |h| & \frac{xprime}{xtilde} \\\\ \n')
for i = 1:length(test_mat)
    t1 = tau{1,2}(i,1); t2 = tau{1,2}(i,2); t3 = absx{1,2}(i,1);
    fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e\\\\ \n',i,...
        t1,t2,t3);
end

fprintf(fid1,'\n');
fprintf(fid1,'\n');
fprintf(fid1,'Terms from error analysis in single precision with diagonal scaling\n');
fprintf(fid1,'Index & \tau & |h| & \frac{xprime}{xtilde} \\\\ \n')
for i = 1:length(test_mat)
    t1 = tau{2,1}(i,1); t2 = tau{2,1}(i,2); t3 = absx{2,1}(i,1);
    fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e\\\\ \n',i,...
        t1,t2,t3);
end

fprintf(fid1,'\n');
fprintf(fid1,'\n');
fprintf(fid1,'Terms from error analysis in double precision without diagonal scaling\n');
fprintf(fid1,'Index & \tau & |h| & \frac{xprime}{xtilde} \\\\ \n')
for i = 1:length(test_mat)
    t1 = tau{2,2}(i,1); t2 = tau{2,2}(i,2); t3 = absx{2,2}(i,1);
    fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e\\\\ \n',i,...
        t1,t2,t3);
end

fclose(fid1);
















