function nrm = norm(X,p)
    if nargin < 2
        p = 2;
    end
    nrm = fp16(norm(double(X),p));
end