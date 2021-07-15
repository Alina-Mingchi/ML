x = [];
for i = 1:20
    x = [x; [randi(10), randi(10)]];
end

for i = 1:20
    x = [x; [randi([50,60]), randi([50,60])]];
end

for i = 1:20
    x = [x; [randi([90,100]), randi([90,100])]];
end

r = randi(3,[200,1]);

%means(run,:) = ones(1, 2); % set up matrix to store means
[k_init, cluster_init, means_init, M_init, distance_init] = Kmeans(3, x, r);
run = 1;
while 1
    [k, cluster, means, M, distance] = Kmeans(3, x, cluster_init);
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

for i = 1:max(cluster)
    plot(x(cluster == i,1),x(cluster == i,2),'o',...
        means(i,1),means(i,2),'*','Color',[rand(),rand(),rand()])
    hold on
end
