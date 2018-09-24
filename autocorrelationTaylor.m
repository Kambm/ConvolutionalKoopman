function [A,dCorrel,derivatives] = autocorrelationTaylor(X,dt,N,nmax)
    % Approximates the autocorrelation of a time series X using
    % a Taylor approximation
    %
    % NOTE : this approximation technique has performed poorly in a number
    % of numerical experiments. It is included for completeness.
    %
    % Parameters
    % ----------
    % X : array-like
    %   n x m matrix of time series data. n is state dimension,
    %   m is time series length.
    % dt : float
    %   time step for time series data
    % N : integer
    %   timesteps of autocorrelation (N*dt = 2*tau)
    % nmax : integer
    %   rank of Taylor approximation
    %
    % Returns
    % -------
    % A : array-like
    %   n*N x n*N matrix of autocorrelation values
    
    [n,m] = size(X);
    A = zeros(n*N,n*N);
    
    % generate derivative filters (precision = 2*order in this case)
    filters{1} = 1;
    for j = 2:nmax+1
        filters{j} = derivativefilter(j-1,dt,j-1);
    end
    
    % generate derivatives
    for j = 1:nmax+1
        derivatives{j} = [];
        for k = 1:n
            nthing = comp(filters{j}', X(k,:));
%             size(nthing)
%             size(derivatives{j})
            derivatives{j} = [derivatives{j}; nthing];
        end
    end
    
    % assemble derivative matrix
    dH = zeros(n*(1+nmax),length(derivatives{nmax}));
    for s = 1:n
        for j = 1:1+nmax
            p = length(derivatives{j})-length(derivatives{nmax});
            dH((s-1)*n+j,:) = derivatives{j}(1+p/2:end-p/2);
        end
    end
    
    %     dH:
    %     [x_1 ]
    %     [x_1']
    %     [... ]
    %     [x_2 ]
    %     [... ]
    dCorrel = zeros(n,n*(nmax+1));
    for j = 0:n-1
        dCorrel(j+1,:) = (dH*dH(N*j+1,:)')';
    end
%     dCorrel = dCorrel/length(derivatives{nmax});
    
    % compute autocovariance
    for j = 1:n
        for k = 1:n
            for p = 1:N
                for q = 1:N
                    for r = 0:2:nmax
                        A(N*(j-1)+p, N*(k-1)+q) = A(N*(j-1)+p,N*(k-1)+q) + dCorrel(j, r+1+(k-1)*(nmax+1))*(dt*(p-q))^r/factorial(r);
                    end
                end
            end
        end
    end
    
    
    
end
