clc;
%%generate list of images 
% add all subdirectories to the Matlab search paths
% addpath(genpath('.'));
% load the pixel data, called "train_pix"
load mfeat-pix.txt -ascii;
%%
% assign data to training and testing sets
n = 1;
m = 1;
for i = 1:2000
    if mod((i/100),2) < 1
        train(n,:) = mfeat_pix(i,:);
        n = n + 1;
    else
        test(m,:) = mfeat_pix(i,:);
        m = m + 1;
    end
end
% %%
m = size(train,1);
n = size(train,2);
 k = 30;

% randomly choose k non-repetitive vector to be the center
center = zeros(k,n);
index = randperm(m,k);
center = train(index,:);
% assign points to clusters, with the criteria of least distance
cluster = zeros(m, 1);
iteration = 0;
Max_iteration = 0;
while true
    % Store old cluster assignments
    old_cluster = cluster;
    % distance: feature vectors
    distance = pdist2(train,center,'euclidean'); 
    Min_distance = zeros(m,1);
    for a = 1:m
        Min_distance(a) = distance(a,1);
        cluster(a) = 1;
    end
    
    for b = 1:m
        for c = 2:k
            if distance(b,c) < Min_distance(b)
                Min_distance(b) = distance(b,c);
                %reassign cluster
                cluster(b) = c;
            end
        end
    end
    
    % Stop when cluster assignments do not change anymore
    if cluster == old_cluster
        break;
    end
    
    % Update the cluster centers
    for d = 1:k
        center(d,:) = mean(train(cluster == d,:));
        
    end
    
    % Stop if number of iteration is more than Max_iteration
    iteration = iteration + 1;
    if iteration > Max_iteration
        break;
    end
end

%%
% visualize centers of clusters
figure('Renderer', 'painters', 'Position', [10 10 900 300])
for g = 1:k
    pic = center(g,:);
    picmatreverse = zeros(15,16);
    picmatreverse(:)= - pic;
    picmat = zeros(15,16);
    for y = 1:15
        subplot(ceil(k/10),10,g);
        picmat(:,y)=picmatreverse(:,16-y);
    end
    pcolor(picmat');
    axis off;
    colormap(gray(10));
end

%%
% create target vectors
for i = 1:10
    for j = 1:100
        class(100*(i-1)+j,:) = [zeros(1, i-1), 1, zeros(1, 10-i)];
    end
end
%%
% calculate weights from linear regression
alpha = 100;
trailing_ones = repmat(1,length(distance),1);
feature_padded = [distance alpha^2*trailing_ones];
class_padded = [class trailing_ones];
C = cov(feature_padded);
I = eye(length(C));

% alternative: without ridge regression
w_opt = (pinv(distance) * class)';

% with ridge regression
weights = (C'*(feature_padded'*class))';
% (feature_padded, target_padded)...
%     .* (C + alpha^2.*I)';
%%
% testing classification outcome for training set
for i = 1:1000
    [M,I] = max(distance(i,:) * w_opt',[],2);
    cat(i,1) = I;
    cat_t(i,:) = distance(i,:) * w_opt';
end
%%
% K-fold cross-validation
K = 10;
validation = zeros(size(train,1)/K,size(train,2));
training = zeros(size(train,1)*(1-1/K),size(train,2));
class_training = zeros(size(train,1)*(1-1/K),10);
class_validation = zeros(size(train,1)/K,10);

% randomly assign data to test vs. validation set
for i = 1:size(train,1)/K
    index = zeros(1, K);
    index = randperm(K); 
    for j = 1:K
       r(j+(i-1)*K) = index(j);  
    end
end


% decide model parameters: 
% model_no = number of models to test
% start = number of minimum features used
model_no = 5;
start = 35;
feature_no = [1,2,5,8,10,15,18,22,25,30,35,40,50,80,100,150,200,240];

% create target vectors
for i = 1:10
    for j = 1:100
        class(100*(i-1)+j,:) = [zeros(1, i-1), 1, zeros(1, 10-i)];
    end
end

MSE_t = zeros(K,model_no);
MSE_v = zeros(K,model_no);
MISS_t = zeros(K,model_no);
MISS_v = zeros(K,model_no);

for j = 1:K
    % designate S_j as validation set
    % designate the union of other S_j' as training set
    
    count_v = 1;
    count_t = 1;
    
    for n = 1:1000
        if r(n) == j
            validation(count_v,:) = train(n,:);
            class_validation(count_v,:) = class(n,:);
            count_v = count_v + 1;
        else
            training(count_t,:) = train(n,:);
            class_training(count_t,:) = class(n,:);
            count_t = count_t + 1;
        end
    end
    
    for it_no = 1:length(feature_no)
        mse_t = 0;
        mse_v = 0;
        miss_t = 0;
        miss_v = 0;
        % extract features for training set
        k = feature_no(it_no);
        distance = zeros(size(training, 1), k);
        m = size(training,1);
        n = size(training,2);
        % randomly choose k non-repetitive vector to be the center
        % center = zeros(k,n);
        index = randperm(m,k);
        center = training(index,:);
        % assign points to clusters, with the criteria of least distance
        cluster = zeros(m, 1);
        iteration = 0;
        Max_iteration = 2;
        while true
            % Store old cluster assignments
            old_cluster = cluster;
            % distance: feature vectors
            distance = pdist2(training,center,'euclidean'); 
            Min_distance = zeros(m,1);
            for a = 1:m
                Min_distance(a) = distance(a,1);
                cluster(a) = 1;
            end
            
            for b = 1:m
                for c = 2:k
                    if distance(b,c) < Min_distance(b)
                        Min_distance(b) = distance(b,c);
                        %reassign cluster
                        cluster(b) = c;
                    end
                end
            end
            
            % Stop when cluster assignments do not change anymore
            if cluster == old_cluster
                break;
            end
            
            % Update the cluster centers
            for d = 1:k
                if ~isempty(find(cluster == d))
                    center(d,:) = mean(training(cluster == d,:));
                   
                end
% %                 center(d,:) = mean(training(cluster == d,:));
% %                 center(~any(~isnan(center),2),:) = [];
            end
            
            
            % Stop if number of iteration is more than Max_iteration
            iteration = iteration + 1;
            if iteration > Max_iteration
                break;
            end
        end
        
%         calculate weights from linear regression
%         alpha = 1;
%         trailing_ones = repmat(1,length(distance),1);
%         feature_padded = [distance alpha^2*trailing_ones];
%         class_padded = [class trailing_ones];
%         C = cov(feature_padded);
%         I = eye(length(C));
        % alternative: without ridge regression
        w_opt = (pinv(distance) * class_training)';
%         with ridge regression
%         weights = (C'*(feature_padded'*class))'
%         (feature_padded, target_padded)...
%             .* (C + alpha^2.*I)';
        
        ind_t = zeros(size(training,1),1);
        % calculate hypothesis matrix for training set
        for i = 1:size(training,1)
            cat_t(i,:) = distance(i,:) * w_opt';
        end
        [M, ind_t] = max(cat_t,[],2);
        [M, ind_class_t] = max(class_training,[],2);
   
        
        % extract features for validation set
        distance_v = zeros(size(validation, 1), k);
        % distance: feature vectors
        distance_v = pdist2(validation,center,'euclidean'); 
        
        % calculate hypothesis matrix for validation set
        for i = 1:size(validation,1)
            cat_v(i,:) = distance_v(i,:) * w_opt';
        end
        [M, ind_v] = max(cat_v,[],2);
        [M, ind_class_v] = max(class_validation,[],2);
        
        
        % MSE, MISS
        for i = 1:size(training,1)
            mse_t = mse_t + immse(cat_t(i), class_training(i));
            if ind_t(i) ~= ind_class_t(i)
                miss_t = miss_t + 1;
            end
        end
        for i = 1:size(validation,1)
            mse_v = mse_v + immse(cat_v(i), class_validation(i));
            if ind_v(i) ~= ind_class_v(i)
                miss_v = miss_v + 1;
            end
        end
        MSE_t(j, it_no) = mse_t/size(training,1);
        MSE_v(j, it_no) = mse_v/size(validation,1);
        MISS_t(j, it_no) = miss_t/size(training,1);
        MISS_v(j, it_no) = miss_v/size(validation,1);
    end
end

MSE_t_mean = mean(MSE_t, 1);
MSE_v_mean = mean(MSE_v, 1);
MISS_t_mean = mean(MISS_t, 1);
MISS_v_mean = mean(MISS_v, 1);

%% plot
figure
plot(feature_no, log10(MSE_t_mean), '--', 'color', 'b')
xlabel('Number of {\it k}-means-based features used')
ylabel('log_{10} mean value from cross-validation run')
hold on
grid on
plot(feature_no, log10(MSE_v_mean), '--', 'color', 'r')
plot(feature_no, log10(MISS_t_mean), 'color', 'b')
plot(feature_no, log10(MISS_v_mean), 'color', 'r')
