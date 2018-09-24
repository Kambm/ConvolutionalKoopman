function A = renormalize(A,S)
    A = inv(S)*A*S;
end