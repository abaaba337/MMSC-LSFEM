function [lev_score, r, time] = leverageScore(A, method, method2)
    % Compute the exact leverage scores of A via rank revealing QR
    arguments
        A
        method = 'Exact' % 'Exact', 'Approx', 'Hashing_Approx'
        method2 = 'svd' % use 'svd' or 'qr'
    end

    [m,n] = size(A);

    if strcmp(method,'Exact')
        if strcmp(method2,'qr')
            tic; [Q,R,~] = qr(A, "econ"); time = toc;
            r = full(sum(abs(diag(R)) > 1e-6)); 
            U = Q(:,1:r);
        elseif strcmp(method2,'svd')
            A = full(A);
            tic; [U,S,~] = svdecon(A); time = toc;
            r = full(sum(diag(S) > 1e-6)); 
            U = U(:,1:r);
        end
        lev_score = full((sum(U.*U, 2))');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif strcmp(method,'Approx')
        A = full(A);
        [S1A, time] = FJLTransform(A);
        S2 = JLTransform(ceil(8*log(m)), n)'; 
        
        if strcmp(method2,'qr')
            tic;
            [~,R,P] = qr(S1A,"econ","vector");
            r = full(sum(abs(diag(R)) > 1e-6)); 
            AP = A(:,P(1:r)); 
            X = R(1:r,1:r)\S2(1:r,:);
            U_approx = AP*X;
            time = time + toc;

        elseif strcmp(method2,'svd')
            tic;
            [~,S,V] = svdecon(S1A); S = diag(S);
            r = full(sum(S>1e-6)); 
            S_inv = sparse(1:r, 1:r, 1./S(1:r));
            R_inv = V(:,1:r)*S_inv;
            X = R_inv*S2(1:r,:);
            U_approx = A*X;
            time = time + toc;
        end
        lev_score = full((sum(U_approx.*U_approx, 2))');

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     elseif strcmp(method,'Hashing_Approx')
   
        s1 = (m > ceil(log(n))) * ceil(log(n)) + (m <= ceil(log(n))) * 2;
        d1 = (m > ceil(n*log(n))) * ceil(n*log(n)) + (m <= ceil(n*log(n))) * ceil((m+n)/2);
        [S1A, time] = sketching(d1, A, 'Hashing', {s1,'Regular'});
        S2 = JLTransform(ceil(8*log(m)), n)';
        
        if strcmp(method2,'qr')
            tic;
            [~,R,P] = qr(S1A, "econ","vector");
            r = full(sum(abs(diag(R)) > 1e-6)); 
            AP = A(:,P(1:r)); 
            R_inv = R(1:r,1:r)\eye(r);
            X = R_inv*S2(1:r,:);
            U_approx = AP*X;
            time = toc + time;

        elseif strcmp(method2,'svd')
            tic;
            S1A = full(S1A);
            [~,S,V] = svdecon(S1A); S = diag(S);
            r = full(sum(S>1e-6)); 
            S_inv = sparse(1:r, 1:r, 1./S(1:r));
            R_inv = V(:,1:r)*S_inv;
            X = R_inv*S2(1:r,:);
            U_approx = A*X;
            time = toc + time;
        end
        lev_score = full((sum(U_approx.*U_approx, 2))');  
    end
end

