function rm = realmax(~)
    e = eps(fp16(1));
    rm = 2^15*(2-e);
end