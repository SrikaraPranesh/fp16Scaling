%%%%%% This script computes the condition number of the scaled matrix and
%%%%%% measures its effectiveness as a preconditioner by computing the
%%%%%% condition number of the preconditioned matrix.

clear all; close all;


rank1_type=1;
addpath('matrices_for_testing')
addpath('GMRES_IR')
addpath('DiagScaling')
addpath('Scaling_Algorithms')
addpath('Misc_Scripts')
rng(1);

%%% Matlab file containing all the test matrices
load test_mat.mat



%%% prec_set=1 is (half,single,double) in GMRES-IR
%%% prec_set=2 is (half,double,quad  ) in GMRES-IR
prec_set=[1;2];

%%% theta -- Headroom to prevent overflow
%%% mu -- Room to prevent underflow
theta=0.1;
mu=1;

%%% scale=11 Performs scaling using Algorithm 2.1
%%% scale=12 Performs scaling using Algorithm 2.2.
%%% scale=2  Performs scaling using just algorithm 2.3.
%%% scale=3 Performs scaling using algorithm 3.2
Scale=3;


%%% flag=1 performs diagonal scaling in Algorithm 3.2
%%% flag=0 no diagonal scaling in Algorithm 3.2
flag=[0;1];

%%% dscale=1, uses Algorithm 2.4 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
%%% dscale=2, uses Algorithm 2.5 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
%%% dscale=3, uses Algorithm 2.6 as diagonal scaling in either Algorithm
%%%           2.3 or Algorithm 3.2
dscale=[1;2;3];


%%%%% condtion number of diagonal scaling algorithms
for prec = 1:2
    for alg = 1:2
        for i=1:length(test_mat)
            
            fprintf('Diagonal scaling part | Matrix %d | Algorithm %d | Precisin %d \n',i,alg,prec);
            
            load(test_mat{i,1});  %%% Load the required matrix
            if (prec_set(prec,1) == 1)
                A=Problem.A;
                A=single(full(A));
                A = single(A);
                n=length(A);
                b=single(randn(n,1));
                u = eps('single');
                [u1,rmins,rmin,rmax,p] = ieee_params('h');
                rmax2 = single(rmax)*single(theta);
            else
                A=Problem.A;
                A=single(full(A));
                n=length(A);
                b=(randn(n,1));
                u = eps('double');
                [u1,rmins,rmin,rmax,p] = ieee_params('h');
                rmax2 = rmax*theta;
            end
            
            Cnumber{prec,alg}(i,1)=cond(A,inf);
            [ A,b,R,C ] = Diagonal_Scaling( A,b,dscale(alg,1),rmax2 );
            Cnumber{prec,alg}(i,2)=cond(A,inf);
            B=double(fp16(A));
            Cnumber{prec,alg}(i,3)=cond(B,inf);
            [L,U,P] = lu(B);
            LL = fp16(double(P')*double(L));
            U=fp16(U);
            At = double((double(U))\((double(L))\((double(P))*(double(A)))));
            Cnumber{prec,alg}(i,4) = cond((double(At)),'inf');
            
        end
        
    end
end


%%%%% condition number of diagonal scaling algorithms
for prec = 1:2
    for alg = 1:2
        
        for i=1:length(test_mat)
            fprintf('rank1 scaling part | Matrix %d | Algorithm %d | Precisin %d \n',i,alg,prec);
            load(test_mat{i,1});  %%% Load the required matrix
            if (prec_set(prec,1) == 1)
                A=Problem.A;
                A=single(full(A));
                A = single(A);
                n=length(A);
                b=single(randn(n,1));
                u = eps('single');
                [u1,rmins,rmin,rmax,p] = ieee_params('h');
                rmax2 = single(rmax)*single(theta);
            else
                A=Problem.A;
                A=single(full(A));
                u = eps('double');
                n=length(A);
                b=randn(n,1);
                [u1,rmins,rmin,rmax,p] = ieee_params('h');
                rmax2 = rmax*theta;
            end
            
            rank1_type=1;
            [ C,b,alpha,beta,C1] = Rank1Update( A,b,theta,flag(alg,1),dscale(1,1),rank1_type);
            Cnumber1{prec,alg}(i,1) = cond(A,inf);
            Cnumber1{prec,alg}(i,2) = cond(C,inf);
            B = double(fp16(C));
            Cnumber1{prec,alg}(i,3) = cond(B,inf);
            
            [L,U,P] = lu(B);
            LL = fp16(double(P')*double(L));
            U = fp16(U);
            At = double((double(U))\((double(L))\( (double(P))*(double(C)))));
            Cnumber1{prec,alg}(i,4) = cond((double(At)),'inf');
            
        end
        
        
    end
end

fid1=fopen('Condition_number.txt','w');



for prec = 1:2
    for alg = 1:2
        
        fprintf(fid1,'Condition numbers diagonal scaling algorithm %d for precision %d \n',alg,prec);
        for i=1:length(test_mat)
            t1 = Cnumber{prec,alg}(i,1); t2 = Cnumber{prec,alg}(i,2);
            t3 = Cnumber{prec,alg}(i,3); t4 = Cnumber{prec,alg}(i,4);
            fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e\\\\ \n',i,...
                t2,t3,t4);
        end
        fprintf(fid1,'\n'); fprintf(fid1,'\n');
    end
end

for prec = 1:2
    for alg = 1:2
        
        fprintf(fid1,'Condition numbers rank1 scaling algorithm %d for precision %d \n',alg,prec);
        for i=1:length(test_mat)
            t1 = Cnumber1{prec,alg}(i,1); t2 = Cnumber1{prec,alg}(i,2);
            t3 = Cnumber1{prec,alg}(i,3); t4 = Cnumber1{prec,alg}(i,4);
            fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e\\\\ \n',i,...
                t2,t3,t4);
        end
        fprintf(fid1,'\n'); fprintf(fid1,'\n');
    end
end

fclose(fid1);
