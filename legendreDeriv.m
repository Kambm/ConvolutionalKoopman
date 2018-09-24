function deriv0 = legendreDeriv(j, n, tau)
    if mod(j-n,2) == 0 && n <= j
        deriv0 = (1/2^j)*sqrt((2*j+1)/2)*(-1)^(0.5*(j-n))*nchoosek(j,0.5*(j-n))*nchoosek(j+n,j)*factorial(n)/(tau^(n+0.5));
    else
        deriv0 = 0;
    end
end