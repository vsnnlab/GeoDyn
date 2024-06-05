clear;
close all;

%% set parameter
% file to open
load('filename.mat');
data_folder = 'data_folder/';

% mask
load('mask.mat'); % load mask image
nanmask = mask;
nanmask(mask==0) = nan;

threshold = 1;
sample_profile = {};
%% Road data
for ff = 1:length(FOI)
    tic;
    load([data_folder, FOI{ff}]);
    % Calculate the Profile
    f_FR = zeros(1,3000);
    f_std = zeros(1,3000);
    for tt = 1:length(f_FR)
        temp = f_sample(:,:,tt).*nanmask;
        f_std(tt) = nanstd(temp(:));
        f_FR(tt) = nanmean(temp(:));
    end

    time_table = auto_sample(f_FR,f_std);
    time_length = time_table(:,2) - time_table(:,1) + 1;
    sample_index = find((time_length>=26) & (time_length<=34));
    time_table = time_table(sample_index,:);
    
    for ss = 1:size(time_table,1)
        count = size(sample_profile,1)+1;
        sample = f_sample(:,:,time_table(ss,1):time_table(ss,2));
        sample(isnan(sample)) = 1e-20;
        sample_profile{count,1} = [data_folder, FOI{ff}];
        sample_profile{count,2} = [time_table(ss,1) time_table(ss,2)];
        temp_temp = imresize(geo_profile(sample,mask,0),[51,30]);
        temp_temp = temp_temp(2:41,:);
        sample_profile{count,3} = temp_temp;
        
        [temp1,temp2] = dyn_profile(sample,mask);
        sample_profile{count,4} = imresize(temp1(:,3:end-2),[37,29]);
        sample_profile{count,5} = imresize(temp2(:,3:end-2),[37,29]);
        
        sample_profile{count,6} = max(sample,[],3).*mask;
        sample_profile{count,7} = latencymap(sample).*mask;
    end
    toc;
end
save filename.mat sample_profile
