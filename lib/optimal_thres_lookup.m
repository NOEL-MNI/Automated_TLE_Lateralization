function [optimal_thres,train_accuracy] = optimal_thres_lookup(tmap,data,group)

    tmax = max(tmap);
    
    for i=1:99

        thres = tmax - tmax*(100-i)/100;
        roi = tmap >= thres;
        
        if ndims(data) == 2
            database = mean(data(:,roi),2);
        else
            database = mean(data(:,roi,:),2);
        end
        
        database = squeeze(database);
        
        [~, accuracy] = repeated_k_fold_lda(database, group, 'diagquadratic',5,24);
        accu(i)=accuracy;
    
    end
    
    [train_accuracy,maxIndex] = max(accu);
    optimal_thres = tmax - tmax*(100-maxIndex)/100;

end
