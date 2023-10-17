function [idx] = iFindZeroCrossing(data)


% data: lines are trials, columns are timeframes

clear val idx idx_02

% This condition virtually never happens. If an exact 0 occurs, replaces it
% by the previous value

for i = 1 : size(data, 1)
    [~, idx_01] = find(data(i,:) == 0);
    data(i,idx_01) = data(i,idx_01 - 1);
end

for i = 1 : size(data,1)
    product_data = data(:,1 : end-1) .* data(:,2 : end);
end

for i = 1 : size(data,1)
    [val idx_02{i}] = find(product_data(i,:) < 0);
end

idx = idx_02;
