function counts = bwlabel_count(labels)   
    label_counts = zeros(max(labels, [], 'all'), 1);  
    
    % First pass
    for (y = 1:size(labels, 1))
        for (x = 1:size(labels, 2))
            if (labels(y,x) > 0)
                label_counts(labels(y,x)) = label_counts(labels(y,x)) + 1;
            end
        end
    end    
    
    % Second pass
    counts = labels;
    for (y = 1:size(counts, 1))
        for (x = 1:size(counts, 2))
            if (counts(y,x) > 0)
                counts(y,x) = label_counts(labels(y,x));
            end
        end
    end
       
end