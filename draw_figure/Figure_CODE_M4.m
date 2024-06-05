clear;
close all;
load('VSDI6_profiles_010518.mat');
load('VSDI_profile_error_010518.mat');
load('Best_ROI_WT_AD_010517.mat')

boundary = 500;
num_data = size(sample_profile,1);
sample_profile_mat_geo = zeros(51,30,num_data);
sample_profile_mat_dyn = zeros(97,29,num_data);

for pp = 1:num_data
    sample_profile_mat_geo(:,:,pp) = imresize(sample_profile{pp,3},[51,30]);
    sample_profile_mat_dyn(:,:,pp) = imresize(sample_profile{pp,5},[97,29]);
end

geo_error_ROI = geo_error(ROI,ROI);
dyn_error_ROI = dyn_error(ROI,ROI);

sample_profile_mat_geo_ROI = sample_profile_mat_geo(:,:,ROI);
sample_profile_mat_dyn_ROI = sample_profile_mat_dyn(:,:,ROI);

cluster_geo = optimal_cluster3(geo_error_ROI,sample_profile_mat_geo_ROI);
cluster_dyn = optimal_cluster3(dyn_error_ROI,sample_profile_mat_dyn_ROI);

max(cluster_geo)
max(cluster_dyn)
%% 4D
temp1 = randperm(500);
temp2 = randperm(500)+500;

train_ROI = sort([temp1(1:400),temp2(1:400)]);
test_ROI = sort([temp1(401:500),temp2(401:500)]);

cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
for ss = 1:boundary*2
    temp_type = (ss>boundary)+1;
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end
cluster_matrix_temp = sum(cluster_matrix_GD,2);
index_temp_1 = cluster_matrix_temp(:,:,1)./(cluster_matrix_temp(:,:,1) + cluster_matrix_temp(:,:,2));
index_temp_1(isnan(index_temp_1)) = 0.5;
index_temp_2 = cluster_matrix_temp(:,:,1) - cluster_matrix_temp(:,:,2);
index_temp = [index_temp_1, index_temp_2];
[~,index_order] = sortrows(index_temp,[1,2]);
cluster_geo = cluster_geo+0.5;
for temp = 1:max(cluster_geo)-0.5
    cluster_geo(cluster_geo==index_order(temp)+0.5) = temp;
end

cluster_matrix_temp = sum(cluster_matrix_GD,1);
index_temp_1 = cluster_matrix_temp(:,:,1)./(cluster_matrix_temp(:,:,1) + cluster_matrix_temp(:,:,2));
index_temp_1(isnan(index_temp_1)) = 0.5;
index_temp_2 = cluster_matrix_temp(:,:,1) - cluster_matrix_temp(:,:,2);
index_temp = [index_temp_1', index_temp_2'];
[~,index_order] = sortrows(index_temp,[1,2]);
cluster_dyn = cluster_dyn+0.5;
for temp = 1:max(cluster_dyn)-0.5
    cluster_dyn(cluster_dyn==index_order(temp)+0.5) = temp;
end

cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
for ss = test_ROI
    temp_type = (ss>boundary)+1;
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end
%%
figure
hold on;
scatter_size = 2;
linnum = 16;
dotspace = 1/(linnum+3);
for gg = 1:max(cluster_geo)
    for dd = 1:max(cluster_dyn)
        temp_index = find(cluster_geo==gg & cluster_dyn==dd);
        if ~isempty(temp_index)
            for ii = 0:size(temp_index)-1
                yy = dd - (floor(ii/linnum)+2)*dotspace;
                xx = gg-1 + (mod(ii,linnum)+2)*dotspace;
                if temp_index(ii+1)<=boundary
                    scatter(xx,yy,scatter_size,'MarkerFaceColor',[0 153 68]/255,'MarkerEdgeColor',[0 153 68]/255);
                else
                    scatter(xx,yy,scatter_size,'MarkerFaceColor',[234 85 20]/255,'MarkerEdgeColor',[234 85 20]/255);
                end
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

%%
scatter_size = 40;
figure
hold on;
scatter(cluster_geo(1:200)-0.75+rand(200,1)*0.5,cluster_dyn(1:200)-0.75+rand(200,1)*0.5,scatter_size,'MarkerFaceColor',[0 153 68]/255,...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
scatter(cluster_geo(201:400)-0.75+rand(200,1)*0.5,cluster_dyn(201:400)-0.75+rand(200,1)*0.5,scatter_size,'MarkerFaceColor',[234 85 20]/255,...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);

for ii = 1:max(cluster_geo)
    plot([ii,ii],[0, max(cluster_dyn)],'k');
end
for ii = 1:max(cluster_dyn)
    plot([0, max(cluster_geo)],[ii,ii],'k');
end
xlim([0 max(cluster_geo)])
ylim([0 max(cluster_dyn)])
axis xy image off

%%
scatter_size = 0.25;
figure
hold on
for gg = 1:max(cluster_geo)
    for dd = 1:max(cluster_dyn)
        sample_size = sum(cluster_matrix_GD(gg,dd,:));
        plot_size = (1/243*sample_size + 242/243)*scatter_size;
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
    plot([ii,ii],[0, max(cluster_dyn)],'k');
end
for ii = 1:max(cluster_dyn)
    plot([0, max(cluster_geo)],[ii,ii],'k');
end
xlim([0 max(cluster_geo)])
ylim([0 max(cluster_dyn)])
axis xy image off



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

%% Find sample
load('mask.mat');
load('rbcmap.mat');
nanmask = mask;
nanmask(nanmask==0) = nan;
gg = 3;
dd = 4;
filename_current = '';
temp_index = find(cluster_geo == gg & cluster_dyn == dd);
mean_profile = mean(sample_profile_mat_dyn_ROI(:,:,temp_index),3);
error_mat = zeros(1,length(temp_index));
for ss = 1:length(temp_index)
    temp = (sample_profile_mat_dyn_ROI(:,:,temp_index(ss)) - mean_profile).^2;
    error_mat(ss) = sum(temp(:));
end
[max_val,ss] = min(error_mat)
error_mat(ss) = inf;
for ss = 11:20
    filename = sample_profile{ROI(temp_index(ss)),1};
    time_point = sample_profile{ROI(temp_index(ss)),2};
    if ~isequal(filename,filename_current)
        load(filename)
    end
    sample = f_sample(:,:,time_point(1):time_point(2));
    time_length = time_point(2)-time_point(1);
    time_bin = time_length/4;
    figure
    for ii = 1:4
        subplot(1,6,ii);
        temp = sample(:,:,round(1+time_bin*(ii-1))).*nanmask;
        hh = imagesc(temp);
        colormap(jet);
        caxis([0 1]);
        axis xy image
        
        set(hh,'alphadata',~isnan(temp))
        set(gca,'color','k');
        set(gcf,'color','w');
        set(gca,'xtick',[],'ytick',[]);
    end
    hh = subplot(1,6,5);
    geo_temp = sample_profile{ROI(temp_index(ss)),3};
    imagesc(geo_temp);
    colormap(hh,hot);
    caxis([0 1])
    axis xy 
    set(gca,'xtick',[],'ytick',[]);
    
    hh = subplot(1,6,6);
    dyn_temp = sample_profile{ROI(temp_index(ss)),4};
    imagesc(dyn_temp);
    colormap(hh,rbcmap);
    caxis([-0.02 0.02])
    colormap(hh,rbcmap);
    axis xy 
    set(gca,'xtick',[],'ytick',[]);
    
    filename_current = filename;
end

%%
figure
for ss = 1:154
    subplot(10,16,ss)
    temp = sample_profile{ROI(temp_index(ss)),4};
    imagesc(temp);
    colormap(rbcmap);
    axis xy off
    title(num2str(ss))
    caxis([-max(abs(temp(:))) max(abs(temp(:)))])
end
