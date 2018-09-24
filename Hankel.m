function H = Hankel(data, maxdelay, delta)
    % Hankel embeds time series data
    %
    % Arguments:
    % - data : n x m time series matrix, where n is the state dimension and
    %   m is the length of the time series.
    % - maxdelay : number of delays incorporated into delay embedded state.
    %   New state has dimension (maxdelay + 1)*n.
    % - delta : spacing between delays.
    %
    % Returns:
    % - H : a generalized Hankel matrix of the time series data. If time
    %   series is one-dimensional, H is a normal Hankel matrix, with
    %   descending rows corresponding to older delays. If time series is
    %   multidimensional, H is a stack of the Hankel matrices for each
    %   component.
    
    s = size(data);
    if s(1) == 1
        H = zeros(maxdelay+1, s(2)-delta*maxdelay);
        for j=0:maxdelay
            H(maxdelay-j+1,:) = data(delta*maxdelay+1-delta*j:length(data)-delta*j);
        end
    else
        H = zeros(s(1)*(maxdelay+1),s(2)-delta*maxdelay);
        for j=1:s(1)
            H(((maxdelay+1)*(j-1)+1):(maxdelay+1)*j,:) = Hankel(data(j,:), maxdelay, delta);
        end
    end
end