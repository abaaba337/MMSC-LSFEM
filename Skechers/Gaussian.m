function [SA, skt_time, t] = Gaussian(d, A, seed, row_per_worker)
    % Gaussian sketching for large A 
    arguments
        d         % d: number of rows for the sketcher
        A         % A: the original system
        seed = [] % random seed
        row_per_worker = [];
    end

    if ~isempty(seed)
        rng(seed)
    end
    
    m = size(A,1); c = 1/sqrt(d);
     
    if ~isempty(seed)
        rng(seed)
    end

    if isempty(row_per_worker)
      row_per_worker = ceil(1e5/m);
    end

    r = mod(d,row_per_worker);
    num_worker = floor(d/row_per_worker);
    new_rows = cell(num_worker,1);
    t = zeros(1,num_worker);
    parfor i = 1:num_worker
        if i == num_worker
            sub_d = row_per_worker + r;
        else
            sub_d = row_per_worker;
        end
        gi = randn(sub_d, m);
        tic;
        new_rows{i} =  c * gi * A ;
        t(i) = toc;
    end
    skt_time = sum(t);
    SA = cell2mat(new_rows);
end
