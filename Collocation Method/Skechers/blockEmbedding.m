function [SB, skt_time, p, I]  = blockEmbedding(d, B, num_worker, method, para)
    arguments
        d          % number of rows for the sketcher
        B          % original system (augmented)
        num_worker % number of parallel workers
        method = 'Gaussian' % type of submatrices embedding
        para = []   % cell of parameters to apply the specific embedding method
                   % 'Hashing'  : {s, hashing_type: 'Regular' or 'Variant'}
                   % 'Levscore' : {'Exact' or 'Approx' or 'Hashing_Approx', 'svd' or 'qr'}

        % p plot(p) to show the visulization for all iterations
        % I is the index for the longest iteration
    end

    [mB, ~] = size(B);
    rowsPerWorker = floor(mB/num_worker);
    rowsPerSubS = max([1,round(d*rowsPerWorker/mB)]);
    subSB = cell(num_worker,1); 
    p = Par(num_worker);

    if strcmp(method,'Gaussian')
        parfor i = 1:num_worker
            if i < num_worker 
                subB = B((i-1)*rowsPerWorker+1:i*rowsPerWorker,:);
            else
                subB = B((i-1)*rowsPerWorker+1:end,:);
            end
            S = sketching_S(rowsPerSubS, subB, 'Gaussian');  
            Par.tic; subSB{i} = S*subB; p(i) = Par.toc;
        end
        stop(p);
    
    elseif strcmp(method,'Hashing')

        if isempty(para)
            para = {1,'Regular'};
        elseif isempty(para{1})
            para{1} = 1;
        elseif isempty(para{2})
            para{2} = 'Regular';
        end
        
        parfor i = 1:num_worker
            if i < num_worker 
                subB = B((i-1)*rowsPerWorker+1:i*rowsPerWorker,:);
            else
                subB = B((i-1)*rowsPerWorker+1:end,:);
            end
            S = sketching_S(rowsPerSubS, subB, 'Hashing', para);
            Par.tic; subSB{i} = S*subB; p(i) = Par.toc;
        end
        stop(p);

    elseif strcmp(method,'Levscore')

        if isempty(para)
            para = {'Exact','svd'};
        elseif isempty(para{1})
            para{1} = 'Exact';
        elseif isempty(para{2})
            para{2} = 'svd';
        end

        if strcmp(para{1},'Exact')
            parfor i = 1:num_worker
                if i < num_worker 
                    subB = B((i-1)*rowsPerWorker+1:i*rowsPerWorker,:);
                else
                    subB = B((i-1)*rowsPerWorker+1:end,:);
                end
                [~, cId] = find(subB(:,1:end-1));
                subCol = subB(:, unique(cId));

                Par.tic;
                [lev_score, ~] = parforleverageScore(subCol, [], [], 'Exact', para{2});
                S = sketching_S(rowsPerSubS, subB, 'Levscore', lev_score/sum(lev_score));
                subSB{i} = S*subB; 
                p(i) = Par.toc;
            end
            stop(p);

        elseif strcmp(para{1},'Approx')
            parfor i = 1:num_worker
                if i < num_worker 
                    subB = B((i-1)*rowsPerWorker+1:i*rowsPerWorker,:);
                else
                    subB = B((i-1)*rowsPerWorker+1:end,:);
                end
                [~, cId] = find(subB(:,1:end-1));
                subCol = subB(:, unique(cId)); [mm,nn] = size(subCol);

                % Construct FJLT, S1 = @(A) S1*A
                D = (rand(mm,1) < 1/2) * 1;  D(D==0) = -1; D = sparse(1:mm, 1:mm, D);
                d = min([ceil(min([(1+(mm/nn))/2,4])*nn), ceil(mm/2)]);
                s = randsample(mm,d);
                S1 = @(A) parforFJLTransform(s, D, A);   

                % Construct S2
                S2 = JLTransform(ceil(8*log(mm)), nn, mm^2); 
                Par.tic;
                [lev_score, ~] = parforleverageScore(subCol, S1, S2', 'Approx', para{2});
                S = sketching_S(rowsPerSubS, subB, 'Levscore', lev_score/sum(lev_score));
                subSB{i} = S*subB; 
                p(i) = Par.toc;
            end
            stop(p);

        elseif strcmp(para{1},'Hashing_Approx')
            parfor i = 1:num_worker
                if i < num_worker 
                    subB = B((i-1)*rowsPerWorker+1:i*rowsPerWorker,:);
                else
                    subB = B((i-1)*rowsPerWorker+1:end,:);
                end
                [~, cId] = find(subB(:,1:end-1));
                subCol = subB(:, unique(cId)); [mm,nn] = size(subCol);

                if mm > ceil(log(nn))
                    s_per_col = ceil(log(nn));
                else
                    s_per_col = 2;
                end
                
                if mm > ceil(log(nn))*nn
                    d = nn*ceil(log(nn)); 
                else
                    d = 2*nn;
                end
                        
                S_1 = sHashing(s_per_col, d, mm, 'Regular'); S1 = @(A) S_1*A; 
                S2  = sHashing(s_per_col, ceil(2*log(mm)), nn, 'Regular');
                Par.tic;
                [lev_score, ~] = parforleverageScore(subCol, S1, S2', 'Approx', para{2});
                S = sketching_S(rowsPerSubS, subB, 'Levscore', lev_score/sum(lev_score));
                subSB{i} = S*subB; 
                p(i) = Par.toc;
            end
            stop(p);
        end

    end
    SB = cell2mat(subSB);
    [skt_time, I] = max([p.ItStop]-[p.ItStart]);
end
