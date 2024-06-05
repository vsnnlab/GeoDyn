function [output_groups,cutoff,cluster_gain,cutoff_list] = optimal_cluster3(input_error,sample_profile_mat)
% sample_profile_mat is 0profile mat*n including n profiles computed in
% input erre

fraggle_lim = 0;
max_cluster_num_temp = inf;
size_lim = max(20,sqrt(size(input_error,1)));
% size_lim = 100000;
while max_cluster_num_temp > size_lim
    Z = linkage(input_error(tril(ones(size(input_error)),-1)==1)','average');
    cutoff_list = Z(:,3);
    cluster_gain = zeros(size(cutoff_list));
    p_0 = mean(sample_profile_mat,3);
    parfor cc = 1:length(cutoff_list)
        cutoff = cutoff_list(cc);
        input_groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
        
        
        %% merge fraggle
        input_cluster = input_groups;
        output_cluster = zeros(size(input_cluster));
        max_cluster = max(input_cluster);
        cluster_order = zeros(max_cluster,3);
        temp = input_cluster;
        
        
        for ii = 1:max_cluster
            cluster_order(ii,1) = mode(temp);
            cluster_order(ii,2) = sum(temp == mode(temp));
            cluster_order(ii,3) = ii;
            temp(temp == mode(temp)) = [];
        end
        others_start = find(cluster_order(:,2)==fraggle_lim,1);
        if isempty(others_start)
            for ii = 1:max_cluster
                output_cluster(input_cluster == cluster_order(ii,1)) = ii;
            end
            max_cluster_num = max_cluster;
        else
            for ii = 1:others_start-1
                output_cluster(input_cluster == cluster_order(ii,1)) = ii;
            end
            
            for ii = others_start:max_cluster
                output_cluster(input_cluster == cluster_order(ii,1)) = others_start;
            end
            max_cluster_num = max_cluster-1;
        end
        
        %% calculate clustering balance
        input_groups = output_cluster;
        max_cluster_num = max(input_groups);
        cluster_gain_temp = zeros(1,max_cluster_num);
        
        for jj = 1:max_cluster_num
            n_j = sum(input_groups==jj);
            p_0_j = mean(sample_profile_mat(:,:,input_groups==jj),3);
            centroid_dist = (p_0-p_0_j).^2;
            centroid_dist = sqrt(sum(centroid_dist(:)));
            cluster_gain_temp(jj) = (n_j-1)*centroid_dist;
        end
        cluster_gain(cc) = sum(cluster_gain_temp);
    end
    
    [~,max_index] = max(cluster_gain);
    cutoff = cutoff_list(max_index);
    output_groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
    
    input_cluster = output_groups;
    output_cluster = zeros(size(input_cluster));
    max_cluster = max(input_cluster);
    cluster_order = zeros(max_cluster,3);
    temp = input_cluster;
    
    for ii = 1:max_cluster
        cluster_order(ii,1) = mode(temp);
        cluster_order(ii,2) = sum(temp == mode(temp));
        cluster_order(ii,3) = ii;
        temp(temp == mode(temp)) = [];
    end
    
    others_start = find(cluster_order(:,2)==fraggle_lim,1);
    if isempty(others_start)
        for ii = 1:max_cluster
            output_cluster(input_cluster == cluster_order(ii,1)) = ii;
        end
        max_cluster_num = max_cluster;
    else
        for ii = 1:others_start-1
            output_cluster(input_cluster == cluster_order(ii,1)) = ii;
        end
        
        for ii = others_start:max_cluster
            output_cluster(input_cluster == cluster_order(ii,1)) = others_start;
        end
        max_cluster_num = max_cluster-1;
    end
    output_groups = output_cluster;
    max_cluster_num_temp = max(output_groups);
    fraggle_lim = fraggle_lim+1;
end