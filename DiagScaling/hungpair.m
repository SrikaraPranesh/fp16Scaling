function [S,u,v] = hungpair(A);
% hungpair computes the Hungarian pair for a complex or real matrix A.
% [S,u,v] = hungpair(A) returns two vectors u,v defining the Hungarian scaling
% H = diag(1./v)*A*diag(1./u) = (1./v)*(1./u)'.*A;
% S is the permutation defining the optimal assignment.
% To put the optimal assignment on the diagonal, use H = H(:,S);
% H is a Hungarian scaled and permuted matrix


%should save on workspace
n = length(A);
B = log(abs(A));
I = -inf*ones(n);
I = diag(ones(n,1))+triu(I,1)+tril(I,-1); %tropical identity

% modify input for algorithm
A(:,:,1,1) = ~isinf(B);
BB = B; BB(isinf(B)) = 0;
A(:,:,1,2) = BB;

A(:,:,2,1) = ~isinf(I);
BB = I; BB(isinf(I)) = 0;
A(:,:,2,2) = BB;

% MAKE MATRIX OF LEAST FINITE TERMS
d = 1;
M = LeastFinite(A,n,d);

% APPLY HUNGARIAN METHOD TO FIND FIRST PERM
[S,SD]=MaxPerXtra(M,n,d);

% now invert permutation
MM = B(:,S);

% now row subtract
M = MM-(diag(MM)*ones(1,n));

% now compute M* by rediculously inefficient method!
MM = M;
MMM = zeros(n); % setup workspace
for k = 1:n
  for i = 1:n
      MMM(i,:) = max(kron(M(i,:)',ones(1,n))+MM,[],1);
  end
  MM=MMM;
end

% U takes its value from max in each col
U = max(MM,[],1);
V = -diag(B(:,S)-kron(ones(n,1),U));
U(S) = U;

u = exp(U); u = u(:);
v = exp(-V); v = v(:);


function M=LeastFinite(A,N,d)

% returns matrix whose components represent the finite monomials of lowest
% degree.

for i=1:1:N
    for j=1:1:N
        M(i,j,1)=0; % finiteness
        M(i,j,2)=0; % degree
        M(i,j,3)=0; % constant
        for k=d+1:-1:1

            if (A(i,j,k,1)==1) % if finite
                M(i,j,1)=1;
                M(i,j,2)=k-1;
                M(i,j,3)=A(i,j,k,2);
            end
        end
    end
end


