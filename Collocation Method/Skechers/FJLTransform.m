function [SHDA, time] = FJLTransform(A, d, close_warning)
    % reference to Achlioptas (2003) and Avron (2010)
    % construct a Fast-Johnson-Lindenstrauss Transform with high probability

    arguments
        A                 % the original system
        d = ceil(min([(1+(size(A,1)/size(A,2)))/2,4])*size(A,2));   % row dim for the new system
        close_warning = true
    end

    [m,n] = size(A); 

    if isempty(d)
        d = ceil(min([(1+(m/n))/2,4])*n);
    end

    if isempty(close_warning)
        close_warning = true;
    end

    d = min([d, round(m/2)]);
    
    d_lb = n*log(m) * log(n*log(m));

    if (d < d_lb) && (~close_warning)
        warning(['FJLTransform: sketching row number too small. ' ...
            'Should have at least have %d rows to construct a %.3f-JLT with high probability.'], ...
            ceil(d_lb), epsilon)
    end

    D = (rand(m,1) < 1/2) * 1; D(D==0) = -1;
    D = sparse(1:m, 1:m, D); DA = full(D*A);
    s = randsample(m,d);
    tic; HDA = dct(DA,'type',4); time = toc;
    SHDA = HDA(s,:);

    % SHDA_col = cell(1,n);
    % parfor j = 1:n
    %     aj = full(A(:,j)); 
    %     HDA = dct(D*aj,'type',4);
    %     SHDA_col{j} = HDA(s,:);
    % end
    % SHDA = cell2mat(SHDA_col);
end

