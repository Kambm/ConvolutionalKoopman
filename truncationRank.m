function r = truncationRank(S)
    % Computes SVD truncation rank by looking for sharpest 'kink'
    % in singular value spectrum. This generally corresponds to the
    % singular values hitting the noise threshold.
    dS = log(diag(S));
    dS = dS(2:end)-dS(1:end-1);
    dS = dS(2:end)-dS(1:end-1);
    [m,r] = max(dS);
end