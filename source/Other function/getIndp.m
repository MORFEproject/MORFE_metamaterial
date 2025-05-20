function Indp = getIndp(n, n_p)
% Find all multi index for n variables and up to n_p order
Indp = cell(n_p,1);

for p_i = 1:n_p
    Indp{p_i} = flip(findAllVectors(n, p_i), 2);
end
end

function R = findAllVectors(n, p)
    R = [];
    R = recursiveFind(0, n, p, [], R);
end

function R = recursiveFind(currentDim, totalDim, remainingSum, currentVector, R)
    if currentDim == totalDim - 1
        R = [R, [currentVector; remainingSum]];
    else
        for i = 0:remainingSum
            R = recursiveFind(currentDim + 1, totalDim, remainingSum - i, [currentVector; i], R);
        end
    end
end
