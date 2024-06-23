function [S_distlevScore, time] = distributedLSS(s,B,num_worker)
    [mB, ~] = size(B);
    rowsPerWorker = floor(mB/num_worker);
    remainder = mod(mB,num_worker);
    S_distlevScore = sparse(s,mB);
    SubS = cell(1,num_worker);
    time = zeros(1,num_worker);
    
    parfor i = 1:num_worker
        if i < num_worker 
            subB = B((i-1)*rowsPerWorker+1:i*rowsPerWorker,:);
        else
            subB = B((i-1)*rowsPerWorker+1:end,:);
        end
       
        rowsPerSubS = max([1,round(s*size(subB,1)/mB)]);
        rowsPerSubS = min([rowsPerWorker+remainder, rowsPerSubS]);
        
        [~, cId] = find(subB);
        subB = subB(:, unique(cId));
        
        tic;
        [lev_score, ~] = leverageScore(subB);
        SubS{i} = leverageScoreSampling(rowsPerSubS, ...
            size(subB,1),lev_score/sum(lev_score));
        time(i) = toc;
    end
    
    r = 1; c = 1;
    for i = 1:num_worker
        [sz1,sz2] = size(SubS{i});
        S_distlevScore(r:r+sz1-1, c:c+sz2-1) = SubS{i};
        r = r + sz1; c = c + sz2;
    end
    S_distlevScore = S_distlevScore/sqrt(num_worker);
end
