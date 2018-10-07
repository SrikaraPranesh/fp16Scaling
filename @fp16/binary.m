function z = binary(y)
% binary(y) is an size(y)-by-18-bit string displaying 
% the s, e, and f fields.
    u = y.u;
    [m,n] = size(u);
    z(m,n) = "";
    for k = 1:m
        for j = 1:n
            v = dec2bin(u(k,j),16);
            z(k,j) = [v(1) ' ' v(2:6) ' ' v(7:16)];
        end
    end
end