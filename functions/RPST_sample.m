function sample_profiles = RPST_sample(sample_profile, Z, cut_off)
cluster_sample_profile = sample_profile;
group_count = ones(size(cluster_sample_profile,1),1);
if cut_off ==1
    index_cutoff = size(Z,1);
else
    index_cutoff = find(Z(:,3)>cut_off,1)-1;
end
for gg = 1:index_cutoff;
    temp_profile_1 = cluster_sample_profile{Z(gg,1)};
    group_count1 = group_count(Z(gg,1));
    cluster_sample_profile{Z(gg,1)} = [];
    temp_profile_2 = cluster_sample_profile{Z(gg,2)};
    group_count2 = group_count(Z(gg,2));
    cluster_sample_profile{Z(gg,2)} = [];
    
    [~,~,scale_profile_1,scale_profile_2] = profile_corr(temp_profile_1,temp_profile_2);
    mean_scale_profile = (scale_profile_1*group_count1 + scale_profile_2*group_count2)/(group_count1+group_count2);
    count = size(cluster_sample_profile,1);
    cluster_sample_profile{count+1} = mean_scale_profile;
    group_count(count+1) = group_count1+group_count2;
end
sample_profiles = cluster_sample_profile;

clear temp;
for ii = 1:size(cluster_sample_profile,1);
    if isempty(cluster_sample_profile{ii})
        temp(ii,1) = 0;
    else
        temp(ii,1) = 1;
    end
end
group_count = group_count.*temp;
temp_group_count = group_count;
plist = find(group_count > 0)';
sample_profiles = cell(length(plist),1);
for pp = 1:length(plist)
    [~,max_index] = max(temp_group_count);
    temp_group_count(max_index) = 0;
    sample_profiles{pp,1} = cluster_sample_profile{max_index};
end