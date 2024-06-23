function S_hashing = sHashing(s,d,m,str)
    if strcmp(str,'Regular')
        r = arrayfun(@(j)randsample(d,s), 1:m, 'UniformOutput', false);
        c = arrayfun(@(j)j*ones(s,1), 1:m, 'UniformOutput', false);
        v = arrayfun(@(j)2.0*binornd(1,0.5,s,1)-1.0, 1:m, 'UniformOutput', false);
        r = cat(1, r{:}); v = cat(1, v{:}); c = cat(1, c{:}); 
        S_hashing = sparse(r, c, v/sqrt(s), d, m);
    
    elseif strcmp(str,'Variant')
        r = arrayfun(@(j)randsample(d,s,true), 1:m, 'UniformOutput', false);
        r_unique = cellfun(@(sample)unique(sample), r, 'UniformOutput', false);  
        c = cellfun(@(rr)ones(length(rr),1), r_unique, 'UniformOutput', false);
        c = arrayfun(@(j)c{j}*j, 1:m, 'UniformOutput', false);
        v = cellfun(@(rr)2.0*binornd(1,0.5,length(rr),1)-1.0, r_unique, 'UniformOutput', false);
        r_unique = cat(1, r_unique{:}); v = cat(1, v{:}); c = cat(1, c{:}); 
        S_hashing = sparse(r_unique, c, v/sqrt(s), d, m);
    end
end
