function [sol, solve_time, sparsity, res] = solveLLS(A, b, sparse_threshold)
    % solve LLS and record time for sparse input A and b 
    arguments 
        A
        b
        sparse_threshold = 0.2 % nnz(A) >= sparse_threshold is considered to be dense
    end

    [m,n] = size(A);
    sparsity = nnz(A)/(m*n);
    if sparsity >= sparse_threshold
        A = full(A);
        tic;
        [Q,R] = qr(A);
        sol = R\(Q'*b);
        solve_time = toc;
    else
        tic;
        [C,R] = qr(A,b);
        sol = R\C;
        solve_time = toc;
    end

    res = norm(A*sol-b);

end

