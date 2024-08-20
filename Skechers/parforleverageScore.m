function [lev_score, r] = parforleverageScore(A, S1, S2, method1, method2)
    % Compute the exact leverage scores of A via rank revealing QR
    arguments
        A  % m*n
        S1 % function to evaluate S1A : S1 = @(A) S1*A
        S2 % S2 : right sketching matrix n * r2 s.t. r2 = O(log(m))
        method1 = 'Exact' % 'Exact' or 'Approx' 
        method2 = 'svd'   %  use 'svd' or 'qr'
    end

    if strcmp(method2,'svd')
        A = full(A);
    end

    if strcmp(method1,'Exact')
        if strcmp(method2,'qr')
            [Q,R,~] = qr(A, "econ");
            r = full(sum(abs(diag(R)) > 1e-6)); 
            U = Q(:,1:r);
        elseif strcmp(method2,'svd')
            [U,S,~] = svdecon(A);
            r = full(sum(diag(S) > 1e-6)); 
            U = U(:,1:r);
        end
        lev_score = full((sum(U.*U, 2))');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif strcmp(method1,'Approx')
        if strcmp(method2,'qr')
            [~,R,P] = qr(S1(A),"econ","vector");
            r = full(sum(abs(diag(R)) > 1e-6)); 
            AP = A(:,P(1:r)); R_inv = R(1:r,1:r)\eye(r);
            X = R_inv*S2(1:r,:);
            U_approx = AP*X;

        elseif strcmp(method2,'svd')
            [~,S,V] = svdecon(S1(A));
            S = diag(S);
            r = full(sum(S>1e-6)); 
            S_inv = sparse(1:r, 1:r, 1./S(1:r));
            R_inv = V(:,1:r)*S_inv;
            X = R_inv*S2(1:r,:);
            U_approx = A*X;
        end
        lev_score = full((sum(U_approx.*U_approx, 2))');
    end
end

