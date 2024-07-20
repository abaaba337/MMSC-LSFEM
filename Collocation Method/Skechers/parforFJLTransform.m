function SHDA = parforFJLTransform(s, D, A)
    % reference to Achlioptas (2003) and Avron (2010)
    % construct a Fast-Johnson-Lindenstrauss Transform with high probability

    arguments
        s                 % id of selected rows
        D                 % random diagonal matrix with independent diagonal entries Dii = +1 with probability 1/2 and Dii = −1 with probability 1/2;
        A                 % the original system
    end

    % [~,n] = size(A); 
    % SHDA_col = cell(1,n);
    % parfor j = 1:n
    %     aj = full(A(:,j)); 
    %     HDA = dct(D*aj,'type',4);
    %     SHDA_col{j} = HDA(s,:);
    % end
    % SHDA = cell2mat(SHDA_col);
    
    DA = full(D*A); 
    HDA = dct(DA,'type',4);
    SHDA = HDA(s,:);
end

