%kmeans clustering, choose random ones as center, calculate the distance
%between each image to the center, and cluster accordingly, then
%recalculate the center, repeating... until the cluster do not change any
%more

%%generate list of images 3
% add all subdirectories to the Matlab search paths
addpath(genpath('.'));
% load the pixel data, resulting in a matlab matrix of dim 2000 x 240
% called "mfeat_pix"
load mfeat-pix.txt -ascii;
mfeat3 = mfeat_pix(601:800,:);
%k = 1;
%k = 2;
k = 3;
%k = 200;
m = size(mfeat3,1);
n = size(mfeat3,2);
center = zeros(k,n);
index = randperm(m,k);
center = mfeat3(index,:);


% The assignments of points to clusters.
cluster = zeros(m, 1);
iteration = 0;
Max_iteration = 10;

while true
    % Store old cluster assignments
    old_cluster = cluster;
    
    % Compute distances from each point to the centers and assign each
    % point to the closest cluster.
    distance = pdist2(mfeat3,center,'euclidean');
    %find minimum distance
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
        center(d,:) = mean(mfeat3(cluster == d,:));
    end
    
    % Stop early if we have performed more than Max_iteration iterations
    iteration = iteration + 1;
    if iteration > Max_iteration
        break;
    end
    
    %visualization
    temp = [mfeat3,cluster];
    vis3 = sortrows(temp,241);
    %if there are too many clusters we visualize the first 10 clusters
    for g = 1:k
        figure(g);
        row_img = find(vis3(:,241) == g);
        num_img = size(row_img);
        visual3 = vis3(:,1:240);
        for h = 1:num_img
            pic = visual3(h,:);
            picmatreverse = zeros(15,16);
            picmatreverse(:)= - pic;
            picmat = zeros(15,16);
            for y = 1:15
                picmat(:,y)=picmatreverse(:,16-y);
            end
            subplot(20,10,h);
            pcolor(picmat');
            axis off;
            colormap(gray(10));
        end
    end    
end
%figure 1 represent all images represented by codebook vector 1



