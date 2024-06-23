function [lev_score, r] = leverageScore(A)
    % Compute the exact leverage scores of A via rank revealing QR
    [Q,R,~] = qr(A, "econ");
    r = full(sum(abs(diag(R)) > 1e-6)); 
    U = Q(:,1:r);
    lev_score = full((sum(U.*U, 2))');
end

