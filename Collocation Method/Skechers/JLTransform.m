function S = JLTransform(d, m, n, epsilon, delta, close_warning)
    % reference to Achlioptas (2003)
    % construct a simple ε-Johnson-Lindenstrauss Transform with probability
    % at least 1 - delta

    arguments
        d   % row dim for the new system
        m   % row dim for the original system
        n   % number of vectors to apply JLT, or col dim for the original system
        epsilon = 0.5
        delta = 0.1
        close_warning = true
    end

    if isempty(epsilon)
        epsilon = 0.5 ;
    end

    d_lb = 6/(epsilon^2) * log(n^2/delta);

    if (d < d_lb) && (~close_warning)
        warning(['JLTransform: sketching row number too small. ' ...
            'Should have at least have %d rows to construct a %.3f-JLT with probability %.3f.'], ...
            ceil(d_lb), epsilon, 1-delta)
    end

    S = sparse(d,m);
    id = find(rand(d,m) < 1/3);
    S(id) = sqrt(3/d);
    id = id(rand(1,length(id)) < 1/2);
    S(id) = -S(id);
end

