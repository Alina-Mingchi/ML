data = importdata('mfeat-pix.txt', ' '); % get data

% d1 and d2: digit range

d1 = 1; % digit (start)
d2 = 200; % digit (end)

% r = randi(k,[200,1]); % random number for assignment
digit = [data(d1:d2,:)]; % get digit data

%% 4-means clustering
k = 4; % number of clusters to obtain

r = randi(k,[200,1]);

%means(run,:) = ones(1, 2); % set up matrix to store means
[k_init, cluster_init, means_init, M_init, distance_init] = Kmeans(k, digit, r);
run = 1;
while 1
    [k, cluster, means, M, distance] = Kmeans(k, digit, cluster_init);
    if cluster == cluster_init
        break
    else
        k_init = k;
        cluster_init = cluster;
        means_init = means;
        M_init = M;
        distance_init = distance;
        run = run + 1;
    end
end

run

figure(1);
for i = 1:k
    pic = means(i,:);  
    picmatreverse = zeros(15,16);
    % the filling of (:) is done columnwise!
    picmatreverse(:)= - pic;
    picmat = zeros(15,16);
    for k = 1:15
        picmat(:,k)=picmatreverse(:,16-k);
    end
    subplot(10,20,i);
    pcolor(picmat');
    axis off;
    colormap(gray(10));
end