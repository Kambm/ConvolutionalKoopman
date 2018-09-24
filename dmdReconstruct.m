function Xdmd = dmdReconstruct(Phi,omega,X0,m,dt)
    % Takes in eigenvectors, eigenvalues and initial condition X0, and
    % performs DMD reconstruction
    %
    % Arguments:
    % - Phi: n x r eigenvector matrix
    % - omega: length-r continuous eigenvalue vector
    % - X0: n x 1 initial condition
    % - dt: timestep
    % - m: number of timesteps in reconstruction
    %
    % Returns:
    % - Xdmd : n x m time series matrix reconstructed from dynamic modes.
    
    b = Phi\X0;
    time_dynamics = zeros(size(Phi,2), m);
    t = (1:m)*dt;
    for iter = 1:m
        time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
    end;
    Xdmd = Phi * time_dynamics;

end