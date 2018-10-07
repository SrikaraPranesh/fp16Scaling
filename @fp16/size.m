function [m,n] = size(X,arg2)
    if nargout == 2
        [m,n] = size(X.u);
    elseif nargin == 2
        m = size(X.u,arg2);
    else
        m = size(X.u);
    end
end