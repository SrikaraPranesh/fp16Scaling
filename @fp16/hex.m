function z = hex(y)
% hex(y) is an 4-bit string.
    u = y.u;
    [m,n] = size(u);
    z(m,n) = "";
    for k = 1:m
        for j = 1:n
            z(k,j) = dec2hex(u(k,j),4);
        end
    end
end