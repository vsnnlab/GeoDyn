function sample_profile = AVG_profile(sample_profile)
profile_size = size(sample_profile{1,1},1);
profile_mat = zeros(profile_size,30,size(sample_profile,1));
for ii = 1:size(sample_profile,1)
    profile_mat(:,:,ii) = imresize(sample_profile{ii,1},[profile_size,30]);
end
sample_profile = mean(profile_mat,3);
end

