function S_levScore = leverageScoreSampling(d,m,p)
    c = zeros(1,d); r = zeros(1,d); v = zeros(1,d);
    for i = 1:d
        r(i) = i;
        [~,c(i)] = find(mnrnd(1, p));
        v(i) = 1/sqrt(d*p(c(i)));
    end
    S_levScore = sparse(r, c, v, d, m);
end

