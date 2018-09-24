function coeff = legendreCoeff(k, j, tau, N)
        coeff = 0;
        for n = 1:N
            if mod(k+n,2) == 0
                innerprod = 0;
            else
                innerprod = 0;
                for l = 0:floor(k/2)
                    innerprod = innerprod + (2^(1-k))*sqrt((2*k+1)/2)*(-1)^l*nchoosek(k,l)*nchoosek(2*k-2*l,k)*tau^(n-0.5)/(k-2*l+n);
                end
            end
            if mod(j-n,2) == 0 && n <= j
                deriv0 = (1/2^j)*sqrt((2*j+1)/2)*(-1)^(0.5*(j-n))*nchoosek(j,0.5*(j-n))*nchoosek(j+n,j)*factorial(n)/(tau^(n+0.5));
            else
                deriv0 = 0;
            end
            coeff = coeff + deriv0*innerprod/factorial(n-1);
        end
end