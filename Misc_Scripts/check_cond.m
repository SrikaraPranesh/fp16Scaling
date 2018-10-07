%CHECK_COND
% function check_cond(A)

% load(['matrices_for_testing\' name])
% A = full(Problem.A);

maxA = norm(abs(A(:)))
rcondA = rcond(A)

A16 = double(fp16(A));
rcondA16 = rcond(A16)
