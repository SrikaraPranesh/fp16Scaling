function plot_scal
[u,rmins,rmin,rmax,p] = ieee_params('h');

tmax = rmax; tmin = rmin;
amin = 1e-10; amax = 1e6;

tmin
fmamin = f(-amin)

fplot(@f,[10*-amin -amin])
grid
% set(gca, 'yscale','log')
% set(gca, 'xscale','log')

function fun = f(x) 
fun = ( (x-amin)*tmax - (x-amax)*tmin ) / (amax - amin);
end

end