function Z = tril(X,arg2)
    if nargin == 1
        arg2 = 0;
    end
    Z = fp16(tril(double(X),arg2));
end