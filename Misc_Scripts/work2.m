%work2    Check permutation invariance.
%Input: n

A = randn(n);
[B,d1,d2] = scale_diag_2side(A)
D1 = diag(d1)
D2 = diag(d2)

P = eye(n); P = P(randperm(n),:);
Q = eye(n); Q = Q(randperm(n),:);

[Bp,d1p,d2p] = scale_diag_2side(P*A*Q)
D1p = diag(d1p)
D2p = diag(d2p)

D1predict = P*D1*P'
D2predict = Q'*D2*Q

err1 = D1p - D1predict
err2 = D2p - D2predict

