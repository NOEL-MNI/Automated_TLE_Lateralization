function [result, accuracy] = k_fold_lda(data, group, type, k)

    partition = cvpartition(group,'KFold',k);
    result = zeros(size(group));
    classification = zeros(size(group));
    
    for i = 1:k
    
        training_data = data(partition.training(i),:);
        training_group = group(partition.training(i));
        
        test_data = data(partition.test(i),:);
        
        model = fitcdiscr(training_data,training_group,'DiscrimType',type);
        classification(partition.test(i)) = predict(model,test_data);
        
    end
    
    result = group == classification;
    accuracy = sum(result)/length(result);

end
