function At = ctranspose(A)
% /ctranspose   change the adjoint flag of A to xor 1

A.adjoint = xor(1, A.adjoint);
At = A;