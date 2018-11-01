% This script computes the number of GMRES iteration and the
% normwise backward error in GMRES based iterative refinement,
% with scaling. For scaling Algorithm 2.3 with Algorithm 2.4 and
% Algorithm 2.5 of the manuscript are used. The iterative refinement
% is performed in two combination of precisions (i) (half,single,double)
% (ii) (half,double,quad). Generates Table 4.3 of the manuscript
%
% Note -- CAUTION!! Read the comments before changing the variables
clear all; close all; clear global variable;
rng(1);

%%% Matlab file containing all the test matrices
load test_mat.mat



%%% prec_set = 1 is (half,single,double) in GMRES-IR
%%% prec_set = 2 is (half,double,quad  ) in GMRES-IR
prec_set = [1;2];


%%% theta -- Headroom to prevent overflow
%%% mu_flag -- Change mu_flag = 1 to 
%%% prevent multiplication by (theta_max X xmax)
%%% in Algorithm 2.3.
theta = 0.1;
mu = 1;
global mu_flag;
mu_flag = 0;

%%% scale = 11 Performs scaling using Algorithm 2.1
%%% scale = 12 Performs scaling using Algorithm 2.2.
%%% scale = 2  Performs scaling using just algorithm 2.3.
Scale = 2 ;

%%% Maximum number of iterative refinement steps
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
            %%% Call the linear system solver
            [x,gmresits{prec,alg}(i,1),irits{prec,alg}(i,1)] = scale64to16solve(A,b,prec_set(prec,1)...
                ,maxit,theta,Scale,dscale(alg,1));
            
            %%% Compute the normwise backward error
            res = double(b) - double(A)*double(x);
            nbe{prec,alg}(i,1) = double(norm(mp(res,34),'inf')/(length(res)*(norm(mp(double(A),34),'inf')*norm(mp(double(x),34),'inf')+ norm(mp(double(b),34),'inf'))));
        end  
    end
end


%%%%%%% PRINT THE LATEX TABLE INTO A FILE %%%%%%%

% creating a text file to print the GMRES iteration table
fid1 = fopen('DiagonalScalingResult.txt','w');
fprintf(fid1,'Column 1 is matrix id \n');
fprintf(fid1,'Columns 2 and 4 are non symmetric diagonal scaling \n');
fprintf(fid1,'Column 3 and 5 are symmetric diagonal scaling \n');
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'Number of GMRES iterations \n');
fprintf(fid1,'Column 2 and 3 are (half,single,double)\n');
fprintf(fid1,'Column 3 and 4 are (half,double,quad)\n');
for i = 1:length(test_mat)
    t1  =  gmresits{1,1}(i,1); t2 = gmresits{1,2}(i,1);
    t11  = irits{1,1}(i,1)-1; t22 = irits{1,2}(i,1)-1;
    t3  =  gmresits{2,1}(i,1); t4 = gmresits{2,2}(i,1);
    t33  =  irits{2,1}(i,1)-1; t44 = irits{2,2}(i,1)-1;
    fprintf(fid1,'%d & %d (%d) & %d (%d) & %d (%d) & %d (%d)\\\\ \n',i,...
        t1,t11,t2,t22,t3,t33,t4,t44);
end

fprintf(fid1,'\n');
fprintf(fid1,'\n');
fprintf(fid1,'backward errors \n');
fprintf(fid1,'Columns 2 and 3 are (half,single,double) \n');
fprintf(fid1,'Columns 4 and 5 are (half,double,quad) \n');
for i = 1:length(test_mat)
    t1 = nbe{1,1}(i,1); t2 = nbe{1,2}(i,1);
    t3 = nbe{2,1}(i,1); t4 = nbe{2,2}(i,1);
    fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e & %6.2e\\\\ \n',i,...
        t1,t2,t3,t4);
end

fclose(fid1);















