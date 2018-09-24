function [x,ya,ta] = legendreModel(xdat,N)
    t = -0.05:0.001:0.05;
    legendreB = @(x,n) legendreP(n,x/0.05)*sqrt(n+1/2)/sqrt(0.05);
    x = [];
    for j = 1:N
        coord = legendreB(t,j-1);
        x(j,:) = comp(legendreB(t,j-1)/norm(coord),xdat);
    end
    x = x';
    A = buildLegendreOperator(N,0.05);
    B = A(1:end-1,end);
    A = A(1:end-1,1:end-1);
    L = 1000:10000;
    sys = ss(A,B,eye(N-1),0*B);
    size(x(L,N));
    size(x(1,1:N-1));
    [ya,ta] = lsim(sys,x(L,N),0.001*(L-1),x(1,1:N-1));
end