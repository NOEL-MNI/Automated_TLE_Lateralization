function [result, accuracy] = repeated_k_fold_lda(data, group, type, k, n)

    result = zeros(length(group),n);
    tmp_accuracy = zeros(n,1);

    parfor i = 1:n
    
        [result(:,i), tmp_accuracy(i)] = k_fold_lda(data, group, type,k);
        
    end
    
    accuracy = mean(tmp_accuracy);

end