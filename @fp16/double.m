function z = double(y)
    u = y.u;
    [m,n] = size(u);
    z = zeros(m,n);
    for k = 1:m
        for j = 1:n
            z(k,j) = unpack16(u(k,j));
        end
    end
    
    % ------------------------------------------------------
    
    function x = unpack16(u)
    % x = unpack16(u) reverses u = pack16(x), x is a double
        sg = bitshift(u,-15);
        s = 1-2*double(sg);   % (-1)^sg
        u = bitxor(u,bitshift(sg,15));
        ebias = bitshift(u,-10);
        u = bitxor(u,bitshift(ebias,10));
        e = double(ebias)-15;
        if e == 16 && u ~= 00 
            f = NaN;
        elseif e == 16
            f = Inf;
        elseif e < -14
            % Denormal
            f = double(u)/1024;
            e = -14;
        else
            % Normal
            f = 1+double(u)/1024;
        end
        x = s*f*2^e;
    end
end
