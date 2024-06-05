clear
close all
load('VSDI_profiles_oldWT4AD3_180409.mat');
load('VSDI_profile_error_oldWT4AD3_180409.mat');
load('Best_ROI_oldWT4AD3_180409.mat')

boundary = 622;
num_data = size(sample_profile,1);
sample_profile_mat_geo = zeros(40,30,num_data);
sample_profile_mat_dyn = zeros(37,29,num_data);
sample_profile_mat_max = zeros(126,97,num_data);
sample_profile_mat_lat = zeros(126,97,num_data);

for pp = 1:num_data
    sample_profile_mat_geo(:,:,pp) = sample_profile{pp,3};
    sample_profile_mat_dyn(:,:,pp) = sample_profile{pp,5};
    sample_profile_mat_max(:,:,pp) = sample_profile{pp,6};
    sample_profile_mat_lat(:,:,pp) = sample_profile{pp,7};
end

geo_error_ROI = geo_error(ROI,ROI)+1e-20;
dyn_error_ROI = dyn_error(ROI,ROI)+1e-20;
max_error_ROI = max_error(ROI,ROI)+1e-20;
lat_error_ROI = lat_error(ROI,ROI)+1e-20;

sample_profile_mat_geo_ROI = sample_profile_mat_geo(:,:,ROI);
sample_profile_mat_dyn_ROI = sample_profile_mat_dyn(:,:,ROI);
sample_profile_mat_max_ROI = sample_profile_mat_max(:,:,ROI);
sample_profile_mat_lat_ROI = sample_profile_mat_lat(:,:,ROI);

%%
figure
WT_geo = squeeze(mean(sample_profile_mat_geo_ROI(:,:,1:length(ROI)/2),3));
imagesc(WT_geo);
axis xy image off
caxis([0 1]);
ylim([1 30]);
colormap(hot);

figure
AD_geo = squeeze(mean(sample_profile_mat_geo_ROI(:,:,length(ROI)/2+1:end),3));
imagesc(AD_geo);
axis xy image off
caxis([0 1]);
ylim([1 30]);
colormap(hot);

figure
imagesc(WT_geo-AD_geo)
axis xy image off
caxis([-0.5 0.5]);
ylim([1 30]);
colormap(jet);
%%
figure
WT_dyn = squeeze(mean(sample_profile_mat_dyn_ROI(:,:,1:length(ROI)/2),3));
imagesc(WT_dyn);
colormap(rbcmap);
caxis([-0.01 0.01]);
axis xy off
figure
AD_dyn = squeeze(mean(sample_profile_mat_dyn_ROI(:,:,length(ROI)/2+1:end),3));
imagesc(AD_dyn);
colormap(rbcmap);
caxis([-0.01 0.01]);
axis xy off
figure
imagesc(WT_dyn-AD_dyn)
axis xy image off
caxis([-0.005 0.005]);
colormap(jet);

%%
geo_sim = geo_error_ROI;
% for ii = 1:size(geo_sim,1)
%     geo_sim(ii,ii) = nan;
% end
% geo_sim = (geo_sim-nanmean(geo_sim(:)))/nanstd(geo_sim(:));
% geo_sim(abs(geo_sim)>2) = nan;
geo_sim = geo_sim-nanmin(geo_sim(:));
geo_sim = geo_sim/nanmax(geo_sim(:));
geo_sim = 1-geo_sim;

temp1 = geo_sim(1:300,1:300);
temp2 = geo_sim(301:600,301:600);
InGroup_geo_sim = [temp1(:); temp2(:)];
InGroup_geo_mean = nanmean(InGroup_geo_sim);
InGroup_geo_std = nanstd(InGroup_geo_sim);

temp = geo_sim(1:300,301:600);
BtGroup_geo_sim = [temp(:);temp(:)];
BtGroup_geo_mean = nanmean(BtGroup_geo_sim);
BtGroup_geo_std = nanstd(BtGroup_geo_sim);

figure
errorbar([InGroup_geo_mean BtGroup_geo_mean],[InGroup_geo_std BtGroup_geo_std],'k.');
xlim([0 3])
ranksum(InGroup_geo_sim, BtGroup_geo_sim)
% ylim([0 1.1])

%%
dyn_sim = dyn_error_ROI;
for ii = 1:size(dyn_sim,1)
    dyn_sim(ii,ii) = nan;
end
dyn_sim = (dyn_sim-nanmean(dyn_sim(:)))/nanstd(dyn_sim(:));
% dyn_sim(abs(dyn_sim)>1) = nan;
dyn_sim = dyn_sim-nanmin(dyn_sim(:));
dyn_sim = dyn_sim/nanmax(dyn_sim(:));
dyn_sim = 1-dyn_sim;

temp = dyn_sim(1:300,1:300);
WT_WT_dyn_mean = nanmean(temp(:));
WT_WT_dyn_std = nanstd(temp(:));

temp = dyn_sim(301:600,301:600);
AD_AD_dyn_mean = nanmean(temp(:));
AD_AD_dyn_std = nanstd(temp(:));

temp = dyn_sim(1:300,301:600);
WT_AD_dyn_mean = nanmean(temp(:));
WT_AD_dyn_std = nanstd(temp(:));

figure
errorbar([WT_WT_dyn_mean AD_AD_dyn_mean WT_AD_dyn_mean],[WT_WT_dyn_std AD_AD_dyn_std WT_AD_dyn_std],'k.');
xlim([0 4])
ylim([0 1.1])
