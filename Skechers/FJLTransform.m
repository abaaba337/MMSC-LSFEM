function [SFDA, time] = FJLTransform(A, d)

    % reference to Achlioptas (2003) and Avron (2010)
    % construct a Fast-Johnson-Lindenstrauss Transform with high probability

    arguments
        A    % the original system
        d = ceil(4*size(A,2)) * (4*size(A,2) <= size(A,1)) ...
          + ceil((size(A,1)+size(A,2))/2) * (4*size(A,2) > size(A,1));    % row dim for the new system
    end

    [m,~] = size(A); 
    D = (rand(m,1) < 1/2) * 1; D(D==0) = -1;
    D = sparse(1:m, 1:m, D); DA = full(D*A);
    s = randsample(m,d);
    tic; FDA = dct(DA,'type',4); time = toc;
    SFDA = FDA(s,:);
end

