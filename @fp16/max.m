function [m,k] = max(x,y)
    if nargin == 1
        [m,k] = max(double(x));
    else
        m = max(double(x),double(y));
    end
    m = fp16(m);
end
