function x = unHankel(H,Ndelay)
    % Extracts the time series data from a Hankel matrix of the
    % data.
    % In the case that H is only approximately Hankel, the
    % reconstructed time series is computed by averaging each
    % time point in the matrix.

    s = size(H);
    if s(1) == Ndelay
        x = zeros(1,s(2)+Ndelay-1);
        for j=1:s(2)+Ndelay-1
            Ah = H(max(j-s(2)+1,1):end, max(1,j-Ndelay+1):min(j,s(2)));
            if size(Ah,2) == 1
                x(j) = Ah(1);
            else
                x(j) = mean(diag(fliplr(Ah)));
            end
        end
    else
        x = zeros(s(1)/Ndelay,s(2)+Ndelay-1);
        for j = 1:s(1)/Ndelay
            x(j,:) = unHankel(H((j-1)*Ndelay+1:Ndelay*j,:), Ndelay);
        end
    end
end