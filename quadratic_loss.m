function L = quadratic_loss(A,B)

L = trace(A\B-eye(size(B,2)))^2;