function [k_new, label, means, M, distance] = Kmeans(k, digit, r)

% k: number of clusters to obtain
% d1 and d2: digit range
% k = 200;
d1 = 1;
d2 = 60;

for j = 1:1:k % identify the cluster assignment for x_i
    xs = 0;
    sum = zeros(1, size(digit,2));
    for i = 1:1:60 % assign clusters for 1st to 200th vector; iterate through all x's
        if r(i,1) == j % if x_i belongs to set S_j
            xs = xs + 1; % number of x's in set S_j
            sum(1,:) = sum(1,:) + digit(i,:); % sum to compute mean for the x's in set S_j
        end
        means(j, :) = sum / xs;
    end
end

distance = ones(size(digit,1),k);
for i = 1:60
    for j = 1:k
        distance(i,j) = norm(digit(i,:) - means(j,:));
    end
end

[M,label] = min(distance,[],2);

for j = 1:k
    if isempty(find(label == j, 1)) == 1
        for i = 1:60
            if label(i,1) > j
                label(i,1) = label(i,1) - 1;
            end
        end
        k = k - 1;
    end
end

k_new = k;
        