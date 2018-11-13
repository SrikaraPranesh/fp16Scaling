% This script computes the number of GMRES iteration and the
% normwise backward error in GMRES based iterative refinement,
% with scaling. For scaling Algorithm 2.1  and
% Algorithm 2.2 of the manuscript are used. The iterative refinement
% is performed in two combination of precisions (i) (half,single,double)
% (ii) (half,double,quad). Generates Table 4.2 of the manuscript
%
% Note -- CAUTION!! Read the comments before changing the variables

clear all;close all;
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
Scale = [11;12] ;

%%% Total number iterative refinement iterations
maxit = 10;

%%% dscale=1, uses Algorithm 2.4 as diagonal scaling in either Algorithm
%%%           2.3 
%%% dscale=2, uses Algorithm 2.5 as diagonal scaling in either Algorithm
%%%           2.3 
dscale = [1;2];

for prec = 1:2
    
    for alg = 1:2
        
        for i = 1:length(test_mat)
            load(test_mat{i,1});  %%% Load the required matrix
            fprintf('I am in matrix number %d \n',i);
            if (prec_set(prec,1)==1)
                A = Problem.A;
                A = single(full(A));
                n = length(A);
                b = single(randn(n,1));
            else
                A = Problem.A;
                A = (full(A));
                n = length(A);
                b = randn(n,1);
            end
            
            [x,gmresits{prec,alg}(i,:),irits{prec,alg}(i,:)] =scale64to16solve(A,b,prec_set(prec,1),...
                                            maxit,theta,Scale(alg,1),dscale);
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
fprintf(fid1,'Number of GMRES iteration \n');
for i = 1:length(test_mat)
    t1 = gmresits{1,1}(i,1); t2 = gmresits{1,2}(i,1);
    t11  = irits{1,1}(i,1)-1; t22 = irits{1,2}(i,1)-1;
    t3 = gmresits{2,1}(i,1); t4 = gmresits{2,2}(i,1);
    t33  =  irits{2,1}(i,1)-1; t44 = irits{2,2}(i,1)-1;
    fprintf(fid1,'%d & %d &(%d) & %d &(%d) & %d &(%d) & %d &(%d)\\\\ \n',i,...
    t1,t11,t2,t22,t3,t33,t4,t44);
end
fclose(fid1);









