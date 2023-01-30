function L = entropy_loss(A,B)

L = trace(A\B)-log(det(A\B))-size(B,2);