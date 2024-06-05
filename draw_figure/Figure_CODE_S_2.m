clear;
close all;
load('VSDI6_profiles_WT_AD_temp_011418.mat');
load('VSDI_profile_error_WT_AD_temp_011418.mat');
% load('Best_ROI_WT_AD_010517.mat')


% ROI = sort(ROI);
Z = linkage(dyn_error,'average');
[H,T,outperm] = dendrogram(Z,0);
% ROI = outperm(1:3815);
ROI = [1:279, 1407:1981];%1:4094; %3815
boundary = 279; %2392; %2113
num_data = size(sample_profile,1);
sample_profile_mat_geo = zeros(51,30,num_data);
sample_profile_mat_dyn = zeros(97,29,num_data);
sample_profile_mat_max = zeros(120,87,num_data);
sample_profile_mat_lat = zeros(120,87,num_data);

for pp = 1:num_data
    sample_profile_mat_geo(:,:,pp) = imresize(sample_profile{pp,3},[51,30]);
    sample_profile_mat_dyn(:,:,pp) = imresize(sample_profile{pp,5},[97,29]);
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

[cluster_geo,cutoff_geo,cluster_gain_geo,cutoff_list_geo] = optimal_cluster3(geo_error_ROI,sample_profile_mat_geo_ROI);
[cluster_dyn,cutoff_dyn,cluster_gain_dyn,cutoff_list_dyn] = optimal_cluster3(dyn_error_ROI,sample_profile_mat_dyn_ROI);
% cluster_max = optimal_cluster3(max_error_ROI,sample_profile_mat_max_ROI);
% cluster_lat = optimal_cluster3(lat_error_ROI,sample_profile_mat_lat_ROI);

max(cluster_geo)
max(cluster_dyn)

% save Rev_fig_4_total_cluster.mat

%% AD WT conditional probability
AD_color = [linspace(255,234,64)',linspace(255,85,64)',linspace(255,20,64)']/255;
WT_color = [linspace(255,0,64)',linspace(255,153,64)',linspace(255,67,64)']/255;
WT_AD_color = cat(1,flipud(WT_color),AD_color);

cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
for ss = 1:length(ROI)
    temp_type = (ROI(ss)>boundary)+1;
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end

%%% AD conditional_prop;
P_AD = (length(ROI)-boundary)/length(ROI);
P_WT = boundary/length(ROI);
P_GD = sum(cluster_matrix_GD,3)'/length(ROI);
temp = cluster_matrix_GD(:,:,2)';
P_GD_AD = temp/sum(temp(:));
temp = cluster_matrix_GD(:,:,1)';
P_GD_WT = temp/sum(temp(:));
P_AD_GD = P_GD_AD.*P_AD./P_GD;
P_WT_GD = P_GD_WT.*P_WT./P_GD;

P_WTAD = (P_AD_GD-P_WT_GD)./(P_AD_GD+P_WT_GD);
P_WTAD(isnan(P_WTAD)) = 0;

%%%
index_temp_1 = nanmean(P_AD_GD,1)';
index_temp_2 = nansum(P_AD_GD,1)';
index_temp = [index_temp_1, index_temp_2];
[~,index_order] = sortrows(index_temp,[1,2],'descend');
cluster_geo = cluster_geo+0.5;
for temp = 1:max(cluster_geo)-0.5
    cluster_geo(cluster_geo==index_order(temp)+0.5) = temp;
end

index_temp_1 = nanmean(P_AD_GD,2);
index_temp_2 = nansum(P_AD_GD,2);
index_temp = [index_temp_1, index_temp_2];
[~,index_order] = sortrows(index_temp,[1,2],'descend');
cluster_dyn = cluster_dyn+0.5;
for temp = 1:max(cluster_dyn)-0.5
    cluster_dyn(cluster_dyn==index_order(temp)+0.5) = temp;
end

cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
for ss = 1:length(ROI)
    temp_type = (ROI(ss)>boundary)+1;
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end

%%%
P_AD = (length(ROI)-boundary)/length(ROI);
P_WT = boundary/length(ROI);
P_GD = sum(cluster_matrix_GD,3)'/length(ROI);
temp = cluster_matrix_GD(:,:,2)';
P_GD_AD = temp/sum(temp(:));
temp = cluster_matrix_GD(:,:,1)';
P_GD_WT = temp/sum(temp(:));
P_AD_GD = P_GD_AD.*P_AD./P_GD;
P_WT_GD = P_GD_WT.*P_WT./P_GD;

P_WTAD = (P_AD_GD-P_WT_GD)./(P_AD_GD+P_WT_GD);
P_WTAD(isnan(P_WTAD)) = 0;


figure
imagesc(P_AD_GD)
colormap(AD_color)
hold on
for ii = 1:max(cluster_geo)
    plot([ii+0.5,ii+0.5],[0+0.5, max(cluster_dyn)+0.5],'color',[0.5 0.5 0.5],'linewidth',0.25,'linestyle',':');
end
for ii = 1:max(cluster_dyn)
    plot([0+0.5, max(cluster_geo)+0.5],[ii+0.5,ii+0.5],'color',[0.5 0.5 0.5],'linewidth',0.25,'linestyle',':');
end
xlim([0 max(cluster_geo)])
ylim([0 max(cluster_dyn)])
axis xy image


figure
imagesc(P_WT_GD)
colormap(WT_color)
hold on
for ii = 1:max(cluster_geo)
    plot([ii+0.5,ii+0.5],[0+0.5, max(cluster_dyn)+0.5],'color',[0.5 0.5 0.5],'linewidth',0.25,'linestyle',':');
end
for ii = 1:max(cluster_dyn)
    plot([0+0.5, max(cluster_geo)+0.5],[ii+0.5,ii+0.5],'color',[0.5 0.5 0.5],'linewidth',0.25,'linestyle',':');
end
xlim([0 max(cluster_geo)])
ylim([0 max(cluster_dyn)])
axis xy image

figure
imagesc(P_WTAD)
colormap(WT_AD_color)
hold on
for ii = 1:max(cluster_geo)
    plot([ii+0.5,ii+0.5],[0+0.5, max(cluster_dyn)+0.5],'color',[0.5 0.5 0.5],'linewidth',0.25,'linestyle',':');
end
for ii = 1:max(cluster_dyn)
    plot([0+0.5, max(cluster_geo)+0.5],[ii+0.5,ii+0.5],'color',[0.5 0.5 0.5],'linewidth',0.25,'linestyle',':');
end
xlim([0 max(cluster_geo)])
ylim([0 max(cluster_dyn)])
axis xy image


%% Pie plot
scatter_size = 0.18;
mult_c = 2.5;
temp = sum(cluster_matrix_GD,3);
num_temp = max(temp(:))-1;
figure
hold on
for gg = 1:max(cluster_geo)
    for dd = 1:max(cluster_dyn)
        sample_size = sum(cluster_matrix_GD(gg,dd,:));
        plot_size = ((mult_c-1)/num_temp*sample_size + (num_temp-mult_c+1)/num_temp)*scatter_size;
        if sample_size ~= 0
            if cluster_matrix_GD(gg,dd,1)>cluster_matrix_GD(gg,dd,1)
                secdraw(gg-0.5,dd-0.5,90,360,plot_size,[0 153 68]/255)
                temp_ratio = cluster_matrix_GD(gg,dd,2)/sample_size;
                if temp_ratio ~= 0
                    secdraw(gg-0.5,dd-0.5,90,360*temp_ratio,plot_size,[234 85 20]/255)
                end
            else
                secdraw(gg-0.5,dd-0.5,90,360,plot_size,[234 85 20]/255)
                temp_ratio = cluster_matrix_GD(gg,dd,1)/sample_size;
                if temp_ratio ~= 0
                    secdraw(gg-0.5,dd-0.5,90,360*temp_ratio,plot_size,[0 153 68]/255)
                end
            end
        end   
    end
end


for ii = 1:max(cluster_geo)
    plot([ii,ii],[0, max(cluster_dyn)],'color',[0.5 0.5 0.5],'linewidth',0.25);
end
for ii = 1:max(cluster_dyn)
    plot([0, max(cluster_geo)],[ii,ii],'color',[0.5 0.5 0.5],'linewidth',0.25);
end
xlim([0 max(cluster_geo)])
ylim([0 max(cluster_dyn)])
axis xy image

%% Individual pie plot

individual_index = [0 279 674 1106 1579 1859 2115 2392 3060 3519 4094];

%
cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),10);
for ss = 1:length(ROI)
    if (ROI(ss)>0 && ROI(ss)<=279)
        temp_type = 1;
    elseif (ROI(ss)>279 && ROI(ss)<=674)
        temp_type = 2;
    elseif (ROI(ss)>674 && ROI(ss)<=1106)
        temp_type = 3;
    elseif (ROI(ss)>1106 && ROI(ss)<=1579)
        temp_type = 4;
    elseif (ROI(ss)>1579 && ROI(ss)<=1859)
        temp_type = 5;
    elseif (ROI(ss)>1859 && ROI(ss)<=2115)
        temp_type = 6;
    elseif (ROI(ss)>2115 && ROI(ss)<=2392)
        temp_type = 7;
    elseif (ROI(ss)>2392 && ROI(ss)<=3060)
        temp_type = 8;
    elseif (ROI(ss)>3060 && ROI(ss)<=3519)
        temp_type = 9;
    elseif (ROI(ss)>3519 && ROI(ss)<=4094)
        temp_type = 10;
    end
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end

discrim_ratio = max(cluster_matrix_GD,[],3)./sum(cluster_matrix_GD,3);
discrim_ratio(isnan(discrim_ratio)) = 0;
figure
imagesc(discrim_ratio)
%%%
scatter_size = 0.2;
mult_c = 2;
temp = sum(cluster_matrix_GD,3);
num_temp = max(temp(:))-1;
color_index = [ 0 255 0; 0 212.5 42.5; 0 170 85; 0 127.5 127.5; 0 85 170; 0 42.5 212.5; 0 0 255; 255 170 0; 255 85 0 ; 255 0 0]/255;
figure
hold on
for gg = 1:max(cluster_geo)
    for dd = 1:max(cluster_dyn)
        sample_size = sum(cluster_matrix_GD(gg,dd,:));
        plot_size = ((mult_c-1)/num_temp*sample_size + (num_temp-mult_c+1)/num_temp)*scatter_size;
        temp_ratio = squeeze(cluster_matrix_GD(gg,dd,:))/sample_size;
        cum_temp_ratio = cumsum(temp_ratio);
        start_angle = [0; cum_temp_ratio(1:end-1)]*360;
        if sample_size ~= 0
            for ii = 10:-1:1
                secdraw(gg-0.5,dd-0.5,0,360*cum_temp_ratio(ii),plot_size,color_index(ii,:))
                
            end
        end
    end
end


for ii = 1:max(cluster_geo)
    plot([ii,ii],[0, max(cluster_dyn)],'k');
end
for ii = 1:max(cluster_dyn)
    plot([0, max(cluster_geo)],[ii,ii],'k');
end
xlim([0 max(cluster_geo)])
ylim([0 max(cluster_dyn)])
axis xy image off


%% Individual probability plot
individual_index = [0 279 674 1106 1579 1859 2115 2392 3060 3519 4094];
color_index = [ 0 255 0; 0 212.5 42.5; 0 170 85; 0 127.5 127.5; 0 85 170; 0 42.5 212.5; 0 0 255; 255 170 0; 255 85 0 ; 255 0 0]/255;
   
%
cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),10);
for ss = 1:length(ROI)
    if (ROI(ss)>0 && ROI(ss)<=279)
        temp_type = 1;
    elseif (ROI(ss)>279 && ROI(ss)<=674)
        temp_type = 2;
    elseif (ROI(ss)>674 && ROI(ss)<=1106)
        temp_type = 3;
    elseif (ROI(ss)>1106 && ROI(ss)<=1579)
        temp_type = 4;
    elseif (ROI(ss)>1579 && ROI(ss)<=1859)
        temp_type = 5;
    elseif (ROI(ss)>1859 && ROI(ss)<=2115)
        temp_type = 6;
    elseif (ROI(ss)>2115 && ROI(ss)<=2392)
        temp_type = 7;
    elseif (ROI(ss)>2392 && ROI(ss)<=3060)
        temp_type = 8;
    elseif (ROI(ss)>3060 && ROI(ss)<=3519)
        temp_type = 9;
    elseif (ROI(ss)>3519 && ROI(ss)<=4094)
        temp_type = 10;
    end
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end


for individual = 1:10
    indi_color = [linspace(1,color_index(individual,1),64)',linspace(1,color_index(individual,2),64)',...
        linspace(1,color_index(individual,3),64)'];
    %%%
%     if individual>6
%         boundary_temp = 0;
%     else
%         boundary_temp = length(temp_index);
%     end
    temp_CM_GD = cluster_matrix_GD(:,:,individual);
    
    P_indi = sum(temp_CM_GD(:))/sum(cluster_matrix_GD(:));
    P_GD = sum(cluster_matrix_GD,3)'/sum(cluster_matrix_GD(:));
    
    P_GD_indi = temp_CM_GD'/sum(temp_CM_GD(:));
    
    P_indi_GD = P_GD_indi.*P_indi./P_GD;
    P_indi_GD(isnan(P_indi_GD)) = 0;

    figure
    imagesc(P_indi_GD)
    colormap(indi_color)
    hold on
    for ii = 1:max(cluster_geo)
        plot([ii+0.5,ii+0.5],[0+0.5, max(cluster_dyn)+0.5],'color',[0.5 0.5 0.5],'linewidth',0.25,'linestyle',':');
    end
    for ii = 1:max(cluster_dyn)
        plot([0+0.5, max(cluster_geo)+0.5],[ii+0.5,ii+0.5],'color',[0.5 0.5 0.5],'linewidth',0.25,'linestyle',':');
    end
    xlim([0 max(cluster_geo)])
    ylim([0 max(cluster_dyn)])
    axis xy image
    caxis([0 max(P_GD_AD(:))])
    colorbar;

    
    %%%
%     scatter_size = 0.15;
%     mult_c = 2;
%     temp = sum(cluster_matrix_GD,3);
%     num_temp = max(temp(:))-1;
% %     num_temp = 243;
%     
%     figure
%     hold on
%     for gg = 1:max(cluster_geo)
%         for dd = 1:max(cluster_dyn)
%             sample_size = sum(cluster_matrix_GD(gg,dd,:));
%             plot_size = ((mult_c-1)/num_temp*sample_size + (num_temp-mult_c+1)/num_temp)*scatter_size;
%             if sample_size ~= 0
%                 if cluster_matrix_GD(gg,dd,1)>cluster_matrix_GD(gg,dd,1)
%                     secdraw(gg-0.5,dd-0.5,90,360,plot_size,color_index(ii,:))
%                     temp_ratio = cluster_matrix_GD(gg,dd,2)/sample_size;
%                     if temp_ratio ~= 0
%                         secdraw(gg-0.5,dd-0.5,90,360*temp_ratio,plot_size,color_index(ii,:))
%                     end
%                 else
%                     secdraw(gg-0.5,dd-0.5,90,360,plot_size,color_index(ii,:))
%                     temp_ratio = cluster_matrix_GD(gg,dd,1)/sample_size;
%                     if temp_ratio ~= 0
%                         secdraw(gg-0.5,dd-0.5,90,360*temp_ratio,plot_size,color_index(ii,:))
%                     end
%                 end
%             end
%         end
%     end
%     
%     
%     for ii = 1:max(cluster_geo)
%         plot([ii,ii],[0, max(cluster_dyn)],'k');
%     end
%     for ii = 1:max(cluster_dyn)
%         plot([0, max(cluster_geo)],[ii,ii],'k');
%     end
%     xlim([0 max(cluster_geo)])
%     ylim([0 max(cluster_dyn)])
%     axis xy image off
end

%% profile axis
figure 
for ii = 1:max(cluster_geo)
    subplot(1,max(cluster_geo),ii)
    temp_index = find(cluster_geo==ii);
%     temp_sample_profile = sample_profile_select(temp_index,3);
    temp_sample_profile = sample_profile_mat_geo_ROI(:,:,temp_index);
    mean_profile = mean(temp_sample_profile,3);
    imagesc(mean_profile);
    colormap(hot);
    
    axis xy image off;
    
    ylim([1 30]);
end

figure
for ii = 1:max(cluster_dyn)
    subplot(max(cluster_dyn),1,ii)
    temp_index = find(cluster_dyn==ii);
    temp_sample_profile = sample_profile_mat_dyn_ROI(:,:,temp_index);
    mean_profile = mean(temp_sample_profile,3);
    imagesc(mean_profile);
    colormap(rbcmap);
    axis xy off;
    caxis([-0.02 0.02])
end

%% Find & plot example sample
% load('mask.mat');
% load('rbcmap.mat');
% nanmask = mask;
% nanmask(nanmask==0) = nan;
% gg = 3;
% dd = 4;
% filename_current = '';
% temp_index = find(cluster_geo == gg & cluster_dyn == dd);
% mean_profile = mean(sample_profile_mat_dyn_ROI(:,:,temp_index),3);
% error_mat = zeros(1,length(temp_index));
% for ss = 1:length(temp_index)
%     temp = (sample_profile_mat_dyn_ROI(:,:,temp_index(ss)) - mean_profile).^2;
%     error_mat(ss) = sum(temp(:));
% end
% [max_val,ss] = min(error_mat);
% error_mat(ss) = inf;
% for ss = 11:20
%     filename = sample_profile{ROI(temp_index(ss)),1};
%     time_point = sample_profile{ROI(temp_index(ss)),2};
%     if ~isequal(filename,filename_current)
%         load(filename)
%     end
%     sample = f_sample(:,:,time_point(1):time_point(2));
%     time_length = time_point(2)-time_point(1);
%     time_bin = time_length/4;
%     figure
%     for ii = 1:4
%         subplot(1,6,ii);
%         temp = sample(:,:,round(1+time_bin*(ii-1))).*nanmask;
%         hh = imagesc(temp);
%         colormap(jet);
%         caxis([0 1]);
%         axis xy image
%         
%         set(hh,'alphadata',~isnan(temp))
%         set(gca,'color','k');
%         set(gcf,'color','w');
%         set(gca,'xtick',[],'ytick',[]);
%     end
%     hh = subplot(1,6,5);
%     geo_temp = sample_profile{ROI(temp_index(ss)),3};
%     imagesc(geo_temp);
%     colormap(hh,hot);
%     caxis([0 1])
%     axis xy 
%     set(gca,'xtick',[],'ytick',[]);
%     
%     hh = subplot(1,6,6);
%     dyn_temp = sample_profile{ROI(temp_index(ss)),4};
%     imagesc(dyn_temp);
%     colormap(hh,rbcmap);
%     caxis([-0.02 0.02])
%     colormap(hh,rbcmap);
%     axis xy 
%     set(gca,'xtick',[],'ytick',[]);
%     
%     filename_current = filename;
% end
% 
% %
% figure
% for ss = 1:154
%     subplot(10,16,ss)
%     temp = sample_profile{ROI(temp_index(ss)),4};
%     imagesc(temp);
%     colormap(rbcmap);
%     axis xy off
%     title(num2str(ss))
%     caxis([-max(abs(temp(:))) max(abs(temp(:)))])
% end


%% GD sort
% cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
% for ss = 1:length(ROI)
%     temp_type = (ROI(ss)>boundary)+1;
%     cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
% end
% cluster_matrix_temp = sum(cluster_matrix_GD,2);
% index_temp_1 = cluster_matrix_temp(:,:,1)./(cluster_matrix_temp(:,:,1) + cluster_matrix_temp(:,:,2));
% index_temp_1(isnan(index_temp_1)) = 0.5;
% index_temp_2 = cluster_matrix_temp(:,:,1) - cluster_matrix_temp(:,:,2);
% index_temp = [index_temp_1, index_temp_2];
% [~,index_order] = sortrows(index_temp,[1,2]);
% cluster_geo = cluster_geo+0.5;
% for temp = 1:max(cluster_geo)-0.5
%     cluster_geo(cluster_geo==index_order(temp)+0.5) = temp;
% end
% 
% cluster_matrix_temp = sum(cluster_matrix_GD,1);
% index_temp_1 = cluster_matrix_temp(:,:,1)./(cluster_matrix_temp(:,:,1) + cluster_matrix_temp(:,:,2));
% index_temp_1(isnan(index_temp_1)) = 0.5;
% index_temp_2 = cluster_matrix_temp(:,:,1) - cluster_matrix_temp(:,:,2);
% index_temp = [index_temp_1', index_temp_2'];
% [~,index_order] = sortrows(index_temp,[1,2]);
% cluster_dyn = cluster_dyn+0.5;
% for temp = 1:max(cluster_dyn)-0.5
%     cluster_dyn(cluster_dyn==index_order(temp)+0.5) = temp;
% end
% 
% cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
% for ss = 1:length(ROI)
%     temp_type = (ROI(ss)>boundary)+1;
%     cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
% end
% 
