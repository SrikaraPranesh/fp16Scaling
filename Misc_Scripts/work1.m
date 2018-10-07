%WORK1

[u,rmins,rmin,rmax,p] = ieee_params('h');
g = u^3;
A = [g g 0;
     0 g g;
     0 0 1/u];