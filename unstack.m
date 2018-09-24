function vs = unstack(v,stackmax)
    statesize = length(v)/stackmax;
    vs = zeros(statesize,stackmax);
    for j = 1:statesize
        for k = 1:stackmax
            vs(j,k) = v((j-1)*stackmax+k);
        end
    end
end
