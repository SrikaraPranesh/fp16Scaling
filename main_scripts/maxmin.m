clear all; close all;
load test_mat.mat

sind = zeros(length(test_mat),1);
for i = 1:length(test_mat)
    load(test_mat{i,1});
    A = Problem.A;
    A = (full(A));
    absA_nozeros = abs(A) + realmax*(A == 0);
    
    rmima(i,1) = length(A);
    rmima(i,3) = max(max(abs(A)));
    rmima(i,4) = min(min(abs(absA_nozeros)));
    rmima(i,2) = cond(A,inf);
    
    if (norm(A-A') == 0)
       sind(i,1)=1; 
    end
    
end