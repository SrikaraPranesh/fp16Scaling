function disp(x)
% fp16 disp
    if isequal(get(0,'format'),'hex')
        disp(hex(x))
    else
        disp(double(x));
    end
end