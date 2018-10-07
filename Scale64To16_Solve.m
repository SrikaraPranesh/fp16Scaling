function [x,gmresits,irits,Cnumber,tau,rx] = Scale64To16_Solve(A,b,sf...
    ,maxit,theta,Scale,flag,dscale)
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
% 3. Make sure that the folder 'matrices_for_testing'
% is in the matlab directory.


if (sf==1)
    A=single(A);
    b=single(b);
    tol=1e-2;
    wp=1;
    rp=2;
    u=eps('single');
    [u1,rmins,rmin,rmax,p] = ieee_params('h');
    rmax2 = single(rmax)*single(theta);
else
    b=b;
    tol=1e-4;
    wp=2;
    rp=4;
    u=eps('double');
    [u1,rmins,rmin,rmax,p] = ieee_params('h');
    rmax2 = rmax*theta;
end

if (Scale ==11 || Scale ==12)
    %%% Performs scaling using Algorithm 2.1 and Algorithm 2.2.
    if (Scale==11)
        A1 = SimpleScaling( A,rmax2,1 );
    elseif (Scale==12)
        A1 = SimpleScaling( A,rmax2,2);
    end
    Cnumber=cond(double(fp16(A1)),inf);
    if (Cnumber <= 1e20)
        [x,gmresits,irits] = gmresir3_1(double(A),double(b),0,wp,rp,maxit,tol,A1);
    else
        x=zeros(length(b),1);
        gmresits=Inf;
        irits=Inf;
    end
    
elseif (Scale==2)
    %%% Performs scaling using just algorithm 2.3.
    Cnumber(1,1)=cond(A,inf);
    A_orig=A;
    b_orig=b;
    [ A,b,R,C ] = Diagonal_Scaling( A,b,dscale,rmax2 );
    Cnumber(1,2)=cond(A,inf);
    Cnumber(1,3)=cond(double(fp16(A)),inf);
    bInf=norm(b,'inf');
    b1=b/bInf;
    if (Cnumber(1,3) <= (1e20) )
%         [x,gmresits,irits]=gmresir3_2(double(A),double(b1),0,wp,rp,maxit,...
%             tol,A_orig,b_orig,C,bInf);
        [x,gmresits,irits]=gmresir3(double(A),double(b1),0,wp,rp,maxit,tol);
        x=bInf*C*x;
    else
        x=zeros(length(b),1);
        gmresits=Inf;
        irits=Inf;
    end
    
elseif (Scale==3)
    %%% Performs scaling using algorithm 3.2, further
    %%% algorithm 2.3 can also be used.
    rank1_type=1;
    [ C,b,alpha,beta,C1] = Rank1Update( A,b,theta,flag,dscale,rank1_type);
    Cnumber(1,1)=cond(A,inf);
    Cnumber(1,2)=cond(C,inf);
    Cnumber(1,3)=cond(double(fp16(C)),inf);
    if (Cnumber(1,3) <= 1e20)
        bInf=norm(b,'inf');
        b1=b/bInf;
        [x,gmresits,irits,tau,rx] = Rank1SystemSolve(C,C1,alpha,beta,b1,wp,rp,maxit,...
                                              tol,sf,bInf);
    else
        x=zeros(length(b),1);
        gmresits(1,1)=Inf;
        gmresits(1,2)=Inf;
        irits(1,1)=Inf;
        irits(1,2)=Inf;
        tau=zeros(1,2);
        rx=0;
    end
    
end

end

