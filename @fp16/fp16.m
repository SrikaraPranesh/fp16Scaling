function y = fp16(x)
% FP16.  Constructor for "fp16" 16-bit floating point,
% also known as "half precision".
% y = fp16(x) has one field, y.u, a uint16 packed with
% one sign bit, 5 exponent bits, and 10 fraction bits.
% See also: fp16/double
% Bug fixes 12/20/2017. See http://blogs.mathworks.com/cleve/2017/12/20.

    if nargin == 0
        y.u = uint16([]);
        y = class(y,'fp16');       
    elseif isa(x,'fp16')
        y = x;
    else
        [m,n] = size(x);
        u = zeros(m,n,'uint16');
        for k = 1:m
            for j = 1:n
                u(k,j) = pack16(double(x(k,j)));
            end
        end
        y.u = u;
        y = class(y,'fp16');
    end

    % ---------------------------------------------------------
    
    function u = pack16(x)
    % u = pack16(x) packs single or double x into uint16 u,
    % with bug fixes 12/20/2017.
        rndevn = @(s) round(s-(rem(s,2)==0.5));
        if x == 0
            u = uint16(0);
        elseif isnan(x)
            u = uint16(Inf);
        elseif isinf(x)
            u = uint16(31744);
        else
            % finite and nonzero x
            [f,e] = log2(abs(x));
            f = 2*f-1;  % Remove hidden bit
            e = e-1;
            if e > 15
                % overflow u, inf = 31*2^10
                u = uint16(31744);
            elseif e < -14
                % denormal u
                u = uint16(rndevn(2^(24+e)*(1+f)));
            else
                % normal
                t = uint16(rndevn(1024*f));
                if t==1024
                    t = uint16(0);
                    e = e+1;
                end
                u = bitxor(t, bitshift(uint16(e+15),10));
            end
        end
        s = 1-min(sign(x)+1,1); % sign bit
        u = bitxor(u,bitshift(s,15));
    end
end
