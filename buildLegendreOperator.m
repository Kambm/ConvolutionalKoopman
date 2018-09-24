function A = buildLegendreOperator(N,tau)
    % Approximates the Koopman operator for on the first N Legendre basis
    % convolutional coordinates. 
    %
    % Arguments:
    % - N : Number of basis vectors.
    % - tau : 1/2 of the delay window length.
    %
    % Returns:
    % - A : an N x N approximation of the Koopman operator on the first N
    %   Legendre convolutional coordinates.
    
    A = zeros(N);
    for k=1:N
        for j=1:N
            A(k,j) = legendreCoeff(k-1,j-1,tau,200);
        end
    end
end