function [A,DB] = buildbasismodel(B,dt,N)
    % Given an orthonormal set of discretized basis functions, approximates
    % the Koopman operator on the convolutional coordinates associated with
    % these basis functions

    % Arguments:
    % - B : an m*N x k matrix of basis vectors, where m is the dimension of
    %   the range of each basis function and N is the number of points of
    %   each function sampled. Each length-N component of these discretized
    %   basis functions should be stacked.
    % - dt : discretization interval (time interval between rows of B).
    % - N : specifies the delay embedding dimension, or equivalently the
    %   number of time points sampled from each basis function.
    
    DB = zeros(size(B));
    for j = 1:size(B,1)/N
        DB((j-1)*N+1:j*N,:) = getDerivatives(B((j-1)*N+1:j*N,:),dt);
    end
    A = B'*DB;
end