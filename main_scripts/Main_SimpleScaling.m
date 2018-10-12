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
Scale = [11;12] ;


%%% flag = 1 performs diagonal scaling in Algorithm 3.2
%%% flag = 0 no diagonal scaling in Algorithm 3.2
flag = 0  ;


%%% Total number iterative refinement iterations
maxit = 10;

%%% dscale = 1, uses Algorithm 2.4 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
%%% dscale = 2, uses Algorithm 2.5 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
%%% dscale = 3, uses Algorithm 2.6 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
dscale = 1;

for prec = 1:2
    
    for alg = 1:2
        
        for i = 1:length(test_mat)
            load(test_mat{i,1});  %%% Load the required matrix
            
            if (prec_set(prec,1)==1)
                A = Problem.A;
                A = single(full(A));
                n = length(A);
                lm(i,1) = n;
                b = single(randn(n,1));
            else
                A = Problem.A;
                A = (full(A));
                n = length(A);
                lm(i,1) = length(A);
                b = (randn(n,1));
            end
            
            [x,gmresits{prec,alg}(i,:),irits{prec,alg}(i,:),Cnumber{prec,alg}(i,:)] =... 
                                            Scale64To16_Solve(A,b,prec_set(prec,1)...
                                            ,maxit,theta,Scale(alg,1),flag,dscale);
                                    
        end
        
        
    end
end


%%%%%%% PRINT THE LATEX TABLE INTO A FILE %%%%%%%

% creating a text file to print the GMRES iteration table
fid1 = fopen('SimpleScalingResult.txt','w');
fprintf(fid1,'Alg 1 is the algorithm where overflow number are mappex to xmax \n');
fprintf(fid1,'Alg 2 is the algorithm where we divide the number by xmax \n');
fprintf(fid1,'Column 1 -- Matrix index \n');
fprintf(fid1,'Column 2 -- Alg 1 (half,single,double)\n');
fprintf(fid1,'Column 3 -- Alg 2 (half,single,double)\n');
fprintf(fid1,'Column 4 -- Alg 1 (half, double, quad)\n');
fprintf(fid1,'Column 5 -- Alg 2 (half, double, quad)\n');
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'Number of GMRES iteration ');
for i = 1:length(test_mat)
    fprintf(fid1,'%d & %d & %d & %d & %d\\\\ \n',i,...
    gmresits{1,1}(i,1),gmresits{1,2}(i,1),gmresits{2,1}(i,1),gmresits{2,2}(i,1));
end
fclose(fid1);









