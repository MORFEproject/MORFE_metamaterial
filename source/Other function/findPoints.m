function selected_indices = findPoints(A, a, error1)
    condition1 = abs(A(2,:) - a) < abs(a) * error1;
    indices1 = find(condition1);
    if ~isempty(indices1)
       indices1 = reduce_indices(indices1',2)';
    end
    selected_indices = A(:,indices1);
end

function reduced_indices = reduce_indices(indices, n)
    indices = sort(indices);
    reduced_indices = [];
    
    start_idx = 1;
    
    for i = 2:length(indices)
        if indices(i) - indices(i - 1) > n
            mid_idx = start_idx + floor((i - start_idx) / 2);
            reduced_indices = [reduced_indices; indices(mid_idx)];
            start_idx = i;
        end
    end
    
    mid_idx = start_idx + floor((length(indices) - start_idx + 1) / 2);
    reduced_indices = [reduced_indices; indices(mid_idx)];
end