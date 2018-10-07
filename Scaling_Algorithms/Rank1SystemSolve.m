function [x,gmresits,irits,tau,rx] = Rank1SystemSolve( C,C1,alpha,beta,b,wp,rp,maxit,tol,sf,...
                                                bInf)
% This algorithm solves linear system of equation
% using scaling of algorithm 3.2 and GMRES based iterative
% refinement.
% C -- Input matrix scaled using algorithm 3.2.
% alpha -- Output of algorithm 3.2.
% beta -- Output of algorithm 3.2.
% b -- Right hand side vector of the linear system.
% A -- Input matrix (Unscaled).
% wp -- Working precision (single or double).
% rp -- Precision in which the residual is computed (double or quad).

n=length(b);
e =(ones(n,1));
u=e;
v=e;
ff=b;

%%%% Solve the system using GMRES based
% %%%% iterative refinement
% [y,gmresits(1,1), irits(1,1)]=gmresir3_2(double(C),double(ff),0,wp,rp,maxit,tol,...
%                               C_org,b_org,C1,bInf);
% [z,gmresits(1,2), irits(1,2)]=gmresir3_2(double(C),double(u),0,wp,rp,maxit,tol,...
%                                C_org,ones(n,1),C1,1);

[y,gmresits(1,1), irits(1,1)]=gmresir3(double(C),double(ff),0,wp,rp,maxit,tol);
[z,gmresits(1,2), irits(1,2)]=gmresir3(double(C),double(u),0,wp,rp,maxit,tol);


if (sf==1)
    y1=double(y);z1=double(z);
else
    y1=y; z1=z;
end

mp.Digits(34);
b11=beta;
f1=mp((mp(b11,34)*mp(v,34)'*mp(y1,34)),34)/mp((1-(mp(b11,34)*mp(e,34)'*mp(z1,34))),34);
f=double(f1);
x = y + (f*z);

%%%% compute quantitites from error analysis
tau(1,1)=double((n*abs(beta)*norm(z1,inf))/mp((1-(mp(b11,34)*mp(e,34)'*mp(z1,34))),34));
tau(1,2)=abs(f);
xp=abs(y1)+((abs(beta*e'*y1)/abs(1-(beta*e'*z)))*abs(z1));
rx=double(norm(xp,inf)/norm(x,inf));

x = alpha*double(x);
x=bInf*(C1)*x;


end

