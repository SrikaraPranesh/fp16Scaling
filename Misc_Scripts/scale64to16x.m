%scale64to16x.m
n = length(A);
theta = 1; mu = 1;

AA = abs(A);
amax = max(AA(:))
amin = min(AA(:)) 

[u,rmins,rmin,rmax,p] = ieee_params('h');

% tmax = min(amax,theta*rmax); tmin = max(amin,mu*rmins);
tmax = theta*rmax;
tmin = mu*rmin;

alpha = (tmax - tmin)/(amax - amin)
beta = (amax*tmin - amin*tmax)/(amax - amin)

C = alpha*A + beta*ones(n)
C16 = fp16(C)
