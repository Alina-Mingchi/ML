clear all;
clc;
%%generate list of images 
% add all subdirectories to the Matlab search paths
addpath(genpath('.'));
% load the pixel data, called "mfeat_pix"
load mfeat-pix.txt -ascii;
%mfeat = mfeat_pix(601:800,:);
%mfeat = mfeat_pix(1:200,:); %This is for digit 0
%mfeat = mfeat_pix(201:400,:); %This is for digit 1
%mfeat = mfeat_pix(401:600,:); %This is for digit 2
%mfeat = mfeat_pix(801:1000,:); %This is for digit 4
%mfeat = mfeat_pix(1001:1200,:); %This is for digit 5
%mfeat = mfeat_pix(1201:1400,:); %This is for digit 6
%mfeat = mfeat_pix(1401:1600,:); %This is for digit 7
%mfeat = mfeat_pix(1601:1800,:); %This is for digit 8
mfeat = mfeat_pix(1801:2000,:); %This is for digit 9
%mfeat = mfeat_pix(1701:1900,:);%Combination of 100 '8'&'9'images respectively
m = size(mfeat,1);
n = size(mfeat,2);
%% K = 1
k = 1;
%randomly choose k non-repetitive vector to be the center
center = zeros(k,n);
index = randperm(m,k);
center = mfeat(index,:);

%assign points to clusters, with the criteria of least distance
cluster = zeros(m, 1);
iteration = 0;
Max_iteration = 10;

while true
    % Store old cluster assignments
    old_cluster = cluster;
    distance = pdist2(mfeat,center,'euclidean');
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
        center(d,:) = mean(mfeat(cluster == d,:));
    end
    
    % Stop if number of iteration is more than Max_iteration
    iteration = iteration + 1;
    if iteration > Max_iteration
        break;
    end
    
    %visualization
    %sort vectors in the same cluster together
    temp = [mfeat,cluster];
    vis3 = sortrows(temp,241);
    for g = 1:k
        figure(g);
        row_img = find(vis3(:,241) == g);
        num_img = size(row_img);
        visual3 = vis3(:,1:240);
        %visualize all vectors in each cluster 
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

%visualize all codebook vectors
figure(k + 1);
for g = 1:k
    pic = center(g,:);
    picmatreverse = zeros(15,16);
    picmatreverse(:)= - pic;
    picmat = zeros(15,16);
    for y = 1:15
        picmat(:,y)=picmatreverse(:,16-y);
    end
    subplot(20,10,g);
    pcolor(picmat');
    axis off;
    colormap(gray(10));
end

%% K = 2
k = 2;
%randomly choose k non-repetitive vector to be the center
center2 = zeros(k,n);
index = randperm(m,k);
center2 = mfeat(index,:);

%assign points to clusters
cluster = zeros(m, 1);
iteration = 0;
Max_iteration = 10;

while true
    % Store old cluster assignments
    old_cluster = cluster;
    distance = pdist2(mfeat,center2,'euclidean');
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
        center(d,:) = mean(mfeat(cluster == d,:));
    end
    
    % Stop if number of iteration is more than Max_iteration
    iteration = iteration + 1;
    if iteration > Max_iteration
        break;
    end
    
    %visualization
    temp = [mfeat,cluster];
    vis3 = sortrows(temp,241);
    m = 0;
    for g = 1:k
        figure(g+10);
        row_img = find(vis3(:,241) == g);
        num_img = size(row_img);
        m = m+num_img(1);
        n = num_img(1);
        visual3 = vis3(:,1:240);
        for h = (m-n+1):m
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

%visualize all codebook vector
figure(k + 11);
for g = 1:k
    pic = center(g,:);
    picmatreverse = zeros(15,16);
    picmatreverse(:)= - pic;
    picmat = zeros(15,16);
    for y = 1:15
        picmat(:,y)=picmatreverse(:,16-y);
    end
    subplot(20,10,g);
    pcolor(picmat');
    axis off;
    colormap(gray(10));
end
%% K = 3
k = 3;
center = zeros(k,n);
index = randperm(m,k);
center = mfeat(index,:);

%assign points to clusters.
cluster = zeros(m, 1);
iteration = 0;
Max_iteration = 10;

while true
    % Store old cluster assignments
    old_cluster = cluster;
    distance = pdist2(mfeat,center,'euclidean');
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
        center(d,:) = mean(mfeat(cluster == d,:));
    end
    
    % Stop if number of iteration is more than Max_iteration
    iteration = iteration + 1;
    if iteration > Max_iteration
        break;
    end
    
    %visualization
    temp = [mfeat,cluster];
    vis3 = sortrows(temp,241);
    c = 0;
    for g = 1:k
        figure(g+20);
        row_img = find(vis3(:,241) == g);
        num_img = size(row_img);
        c = c+num_img(1);
        d = num_img(1);
        visual3 = vis3(:,1:240);
        for h = (c-d+1):c
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

%visualize all codebook vectors
figure(k + 21);
for g = 1:k
    pic = center(g,:);
    picmatreverse = zeros(15,16);
    picmatreverse(:)= - pic;
    picmat = zeros(15,16);
    for y = 1:15
        picmat(:,y)=picmatreverse(:,16-y);
    end
    subplot(20,10,g);
    pcolor(picmat');
    axis off;
    colormap(gray(10));
end

%% K = 200
k = 200;
center = zeros(k,n);
index = randperm(m,k);
center = mfeat(index,:);

%assign points to clusters.
cluster = zeros(m, 1);
iteration = 0;
Max_iteration = 10;

while true
    % Store old cluster assignments
    old_cluster = cluster;
    distance = pdist2(mfeat,center,'euclidean');
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
    
    %Update the cluster centers
    new_center = zeros(200,240);
    for d = 1:k      
        new_center(d,:) = mean(mfeat(cluster == d,:));
    end
    
    % Stop if number of iteration is more than Max_iteration
    iteration = iteration + 1;
    if iteration > Max_iteration
        break;
    end
    
    %visualization
    temp = [mfeat,cluster];
    vis3 = sortrows(temp,241);
    %we plot the first 10 clusters
    for g = 1:10
        figure(g+30);
        row_img = find(vis3(:,241) == g);
        visual3 = vis3(:,1:240);
        pic = visual3(g,:);
        picmatreverse = zeros(15,16);
        picmatreverse(:)= - pic;
        picmat = zeros(15,16);
        for y = 1:15
            picmat(:,y)=picmatreverse(:,16-y);
        end
        subplot(20,10,g);
        pcolor(picmat');
        axis off;
        colormap(gray(10));
    end    
end

%visualize first 10 codebook vectors
figure(k + 31);
for g = 1:10
    pic = center(g,:);
    picmatreverse = zeros(15,16);
    picmatreverse(:)= - pic;
    picmat = zeros(15,16);
    for y = 1:15
        picmat(:,y)=picmatreverse(:,16-y);
    end
    subplot(20,10,g);
    pcolor(picmat');
    axis off;
    colormap(gray(10));
end
