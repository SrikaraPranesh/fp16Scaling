function D = diag(X,k)
    if nargin == 1
        k = 0;
    end
    D = fp16(diag(double(X),k));
end