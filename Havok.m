function [H, U, S, V, x, dx, sys, y, t] = Havok(xdat,stackmax,dt,rmax)
    % Applies the full stack HAVOK procedure to a time series xdat.
    %
    % Arguments:
    % - xdat : time series data. Must be sampled uniformly in time. Columns
    %   correspond to state snapshots.
    % - stackmax : maximum delay number.
    % - dt : time interval between snapshots in xdat.
    % - rmax : maximum allowed truncation rank.
    %
    % Returns:
    % - H : Hankel matrix of the time series xdat.
    % - U : matrix of right singular vectors of H.
    % - S : diagonal matrix of singular values of H.
    % - V : matrix left singular vectors of H.
    % - x : matrix of convolutional coordinate history.
    % - dx : matrix of time derivative of convolutional coordinates.
    % - sys : linear system with forcing coordinate identified using HAVOK
    %   procedure.
    % - y : approximation of x using HAVOK model
    % - t : time range of HAVOK approximation
    
    
    clear('U','S','V','H');
    addpath('./utils');
    % COMPUTE
    H = Hankel(xdat,stackmax,1);

    [U,S,V] = svd(H,'econ');

    % Obtain truncation rank via inspection of singular value spectrum.
    % Selects rank with maximal second derivative, generally corresponding
    % to where the singular values reach machine precision or the noise
    % threshold. Tune as needed.
    logSigs = log(diag(S));
    logSigs = logSigs(2:end)-logSigs(1:end-1);
    logSigs = logSigs(2:end)-logSigs(1:end-1);
    [m,r] = max(logSigs);
    r=min(rmax,r);
    
    V = (S*V')';

    
    % compute derivative using fourth order central difference
    % use TVRegDiff if more error 
    dV = zeros(length(V)-5,r);
    for i=3:length(V)-3
        for k=1:r
            dV(i-2,k) = (1/(12*dt))*(-V(i+2,k)+8*V(i+1,k)-8*V(i-1,k)+V(i-2,k));
        end
    end  
    % concatenate
    x = V(3:end-3,1:r);
    dx = dV;
    
    % Compute least squares regression model
    A = x\dx;
    % Separate model and forcing term
    B = A(:,r);
    A = A(:,1:r-1);
    
    L = 1:length(x);
    sys = ss(A,B,eye(r-1),0*B);
    size(A)
    size(B)
    size(x)
    [y,t] = lsim(sys,x(L,r),dt*(L-1),x(1,1:r-1));
end