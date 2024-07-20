function S_hashing = sHashing(s, d, m, str)
    
    % d : sketched row number
    % m : original row number
    % s : s non-zeros per column
    
    arguments
        s
        d
        m
        str = 'Regular'
    end
    
    if strcmp(str,'Regular')
        r = cell(m, 1);
        parfor j = 1:m
            r{j} = randsample(d,s);
        end
        r = cell2mat(r);
        c = kron((1:m)',ones(s,1));
        v = 2*binornd(1,0.5,s*m,1)-1;
        S_hashing = sparse(r, c, v/sqrt(s), d, m);
    
    elseif strcmp(str,'Variant')
        r = cell(m, 1);
        c = cell(m, 1);
        parfor j = 1:m
            rj = unique(randsample(d,s,true));
            r{j} = rj;
            c{j} = ones(length(rj),1)*j;
        end
        r = cell2mat(r);
        c = cell2mat(c);
        v = 2*binornd(1,0.5,length(c),1)-1;
        S_hashing = sparse(r, c, v/sqrt(s), d, m);
    end
end
