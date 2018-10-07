function z = eps(x)
   [~,e] = log2(abs(double(x)));
   z = fp16(pow2(1,e-11));
end