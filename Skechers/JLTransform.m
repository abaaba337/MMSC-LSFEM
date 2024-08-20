function S = JLTransform(d, m)
    % reference to Achlioptas (2003)
    % construct a simple Johnson-Lindenstrauss Transform 

    arguments
        d   % row dim for the new system
        m   % row dim for the original system
    end

    S = sparse(d,m);
    id = find(rand(d,m) < 1/3);
    S(id) = sqrt(3/d);
    id = id(rand(1,length(id)) < 1/2);
    S(id) = -S(id);
end

