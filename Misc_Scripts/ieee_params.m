function [u,rmins,rmin,rmax,p] = ieee_params(prec)
% IEEE_PARAMS   Parameters for IEEE arithmetic.
%   [u,rmins,rmin,rmax,p] = IEEE_PARAMS(prec) returns 
%      u: the unit roundoff,
%      rmins: the smallest positive floating-point number (subnormal),
%      rmin: the smallest positive normalized floating-point number,
%      rmax: the largest floating-point number,
%      p: the number of binary digits in the significand,
%    where prec = 'half', 'single', 'double' (the default), 'quadruple'.
%    With no input and output arguments, IEEE_PARAMS prints a table showing
%    all the parameters for all four precisions.
%    Note: the max and min are not representable in double precison for 'quad'.

% Reference:
% IEEE Standard for Floating-Point Arithmetic, IEEE Std 754-2008 (revision 
% of IEEE Std 754-1985), 58, IEEE Computer Society, 2008; pages 8, 13.  

if nargin < 1, prec = 'd'; end

switch prec(1)
  case 'h'
    p = 11; emax = 15;
  case 's'
    p = 24; emax = 127;
  case 'd'
    p = 53; emax = 1023;
  case 'q' % Quadruple included but max and min not representable!
    p = 113; emax = 16383;
end
    
emin = 1-emax; % For all formats.
rmins = 2^emin * 2^(1-p);
rmin = 2^emin;
rmax = 2^emax * (2-2^(1-p));
u = 2^(-p);

if nargin < 1 && nargout < 1
   precs = 'hsdq';
   fprintf('        u        rmins       rmin       rmax    p\n')
   fprintf('    ------------------------------------------------\n')
   for j = 1:4
      [u,rmins,rmin,rmax,p] = ieee_params(precs(j));
      fprintf('%s: %9.2e  %9.2e  %9.2e  %9.2e  %3.0f\n',...
               precs(j),u,rmins,rmin,rmax,p)
   end
end   
