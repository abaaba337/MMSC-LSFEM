function [lev_score, r, time] = leverageScore(A, method, para, method2)
    % Compute the exact leverage scores of A via rank revealing QR
    arguments
        A
        method = 'Exact' % 'Exact', 'Approx', 'Hashing_Approx'
        para = {[], [], []} 
        % para for FJLTransform if Approx: {d, close_warning}
        % para for sHashing if Hashing_Approx: {s,d,str}
        method2 = 'svd' % use 'svd' or 'qr'
    end

    [m,n] = size(A);

    if isempty(para)
        para = {[], [], []};
    end

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
        [S1A, time] = FJLTransform(A, para{1}, para{2});
        S2 = JLTransform(ceil(8*log(m)), n, m^2)'; 
        
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

         if isempty(para{1})
             if m > ceil(log(n))
                 para{1} = ceil(log(n));
             else
                 para{1} = 2;
             end
         end

         if isempty(para{2})
             if m/n > ceil(log(n))
                 para{2} = n*ceil(log(n));
             else
                 para{2} = round(2*n);
             end
         end

        if isempty(para{3})
             para{3} = 'Regular';
        end

        [S1A, skt_time1] = sketching(para{2}, A, 'Hashing', para([1,3]));
        
        if strcmp(method2,'qr')
            tic;
            [~,R,P] = qr(S1A, "econ","vector");
            r = full(sum(abs(diag(R)) > 1e-6)); 
            AP = A(:,P(1:r)); 
            R_inv = R(1:r,1:r)\eye(r);
            time = toc + skt_time1;
            [Xt, skt_time2] = sketching(ceil(2*log(m)), R_inv', 'Hashing', ...
                para([1,3])); % 2*r
            tic; U_approx = AP*Xt';  time = toc + time + skt_time2;
        elseif strcmp(method2,'svd')
            tic;
            S1A = full(S1A);
            [~,S,V] = svdecon(S1A); S = diag(S);
            r = full(sum(S>1e-6)); 
            S_inv = sparse(1:r, 1:r, 1./S(1:r));
            R_inv = V(:,1:r)*S_inv;
            time = toc + skt_time1;
            [Xt, skt_time2] = sketching(ceil(2*log(m)), R_inv', 'Hashing', ...
                para([1,3]));
            tic; U_approx = A*Xt';  time = toc + time + skt_time2;
        end
        lev_score = full((sum(U_approx.*U_approx, 2))');  
    end
end

