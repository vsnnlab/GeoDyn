function [groups,Z,cutoff,CVI_temp,cutoff_list] = optimal_cluster2(input_error)
temp = ones(size(input_error));
dissim = input_error(tril(temp,-1)==1)';
Z = linkage(dissim,'average');
cutoff_list = Z(:,3);
CVI_temp = zeros(size(cutoff_list));
Intra_Clust_Energy = zeros(size(cutoff_list));
Inter_Clust_Energy = zeros(size(cutoff_list));

parfor cc = 1:length(cutoff_list)
    cutoff = cutoff_list(cc);
    groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
%     CVI_temp(ii) = Cluster_Balance(temp_groups,input_error);
    input_groups = groups;
    
    max_cluster_num = max(groups);

    %% calculate clustering balance
    IntraCE_ii = zeros(1,max_cluster_num);

    for ii = 1:max_cluster_num
        temp_corr = input_error(input_groups == ii,input_groups == ii);
        dissim_temp_corr = temp_corr(find(tril(ones(size(temp_corr)))==0));
        IntraCE_ii(ii) = sum(dissim_temp_corr.^2)/sum(input_groups == ii);
    end
    Intra_Clust_Energy(cc) = sum(IntraCE_ii);
    
    %inter distance
    InterCE_ij = zeros(max_cluster_num,max_cluster_num);
    for ii = 1:max_cluster_num
        for jj = 1:max_cluster_num
            temp_corr_jj = input_error(input_groups == ii,input_groups == jj);
            dissim_tri_temp_corr_jj = temp_corr_jj(find(tril(ones(size(temp_corr_jj)))==1));
            InterCE_ij(ii,jj) = sum(dissim_tri_temp_corr_jj(:))/sum(input_groups == ii)/sum(input_groups == jj);
        end
    end
    tri_InterCE_ij = InterCE_ij(logical(tril(ones(size(InterCE_ij)),-1)));
    Inter_Clust_Energy(cc) = sum(tri_InterCE_ij.^2)/max(input_groups);
    
    CVI_temp(cc) = Intra_Clust_Energy(cc)+Inter_Clust_Energy(cc);
    
end
 CVI_temp= Intra_Clust_Energy+Inter_Clust_Energy;

[~,min_index] = min(CVI_temp);
cutoff = cutoff_list(min_index);
groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
[~,~,outperm] = dendrogram(Z,0,'colorthreshold',cutoff);
    % group renumbering
    groups = groups+0.5;
    temp_groups = groups(outperm);
    new_num = 1;
    while ~isempty(temp_groups)
        target_num = temp_groups(1);
        groups(groups == target_num) = new_num;
        temp_groups(temp_groups==target_num) = [];
        new_num = new_num+1;
    end