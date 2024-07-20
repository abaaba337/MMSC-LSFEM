function S_levScore = leverageScoreSampling(d,m,p)
    r = 1:d; 
    c = randsample(m, d, true, p);
    v = 1./sqrt(d.*p(c));
    S_levScore = sparse(r, c, v, d, m);
end

