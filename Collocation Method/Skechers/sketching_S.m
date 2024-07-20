function S = sketching_S(d, A, method, para)

    arguments
        d         % d: number of rows for the sketcher
        A         % A: the original system
        method    % 'Gaussian', 'Hashing', 'Levscore'
        para = []      

        % 'Gaussian' para : seed
        % 'Hashing'  para : {s, hashing_type: 'Regular' or 'Variant'}
        % 'Levscore' para : p
    end

    m = size(A,1);

    if strcmp(method,'Gaussian')

        if ~isempty(para)
            rng(para)
        end

        S = randn(d, m)/sqrt(d);

    elseif strcmp(method,'Hashing')

        if isempty(para)
            para = {1,'Regular'};
        elseif isempty(para{1})
            para{1} = 1;
        else 
            para{2} = 'Regular';
        end

        S = sHashing(para{1}, d, m, para{2});

    elseif strcmp(method,'Levscore')
        S = leverageScoreSampling(d,m,para);
    end

end

