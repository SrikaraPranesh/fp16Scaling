%WORK4
%Input: alpha

A = [1 1 alpha
     1 -1 alpha
     1  1 0]
[B,d1prod,d2prod,its] = scale_diag_2side_symm(A)
