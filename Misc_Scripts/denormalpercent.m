function [ dnp ] = denormalpercent(B,rmin,rmins )
%DENORMALPERCENT percentage of denormal numbers in the matrix
%   This function check the what percentage of numbers become
%   subnormal when converted to fp16. 


[x1,y1] = find(B < rmin & B >= rmins);

tn = length(B)*length(B);

dnp = (abs(length(x1))/(tn))*100;

end

