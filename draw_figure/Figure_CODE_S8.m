%%
temp = randperm(4149);
train_ROI = temp(1:2000);
test_ROI = temp(2001:4000);

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


cluster_matrix_M = zeros(max(cluster_max),2);
cluster_matrix_L = zeros(max(cluster_lat),2);

for ss = 1:length(train_ROI)
    temp_type = (train_ROI(ss)>2113)+1;
    cluster_matrix_M(cluster_max(train_ROI(ss)), temp_type) =  cluster_matrix_M(cluster_max(train_ROI(ss)), temp_type)+1;
    cluster_matrix_L(cluster_lat(train_ROI(ss)), temp_type) =  cluster_matrix_L(cluster_lat(train_ROI(ss)), temp_type)+1;
end

% sorting
cluster_matrix_temp = cluster_matrix_M;

index_temp_1 = cluster_matrix_temp(:,1)./(cluster_matrix_temp(:,1) + cluster_matrix_temp(:,2));
index_temp_1(isnan(index_temp_1)) = 0.5;
index_temp_2 = cluster_matrix_temp(:,1) - cluster_matrix_temp(:,2);
index_temp = [index_temp_1, index_temp_2];

[~,index_order] = sortrows(index_temp,[1,2]);

cluster_max = cluster_max+0.5;
for draw_deck = 1:max(cluster_max)-0.5
    cluster_max(cluster_max==index_order(draw_deck)+0.5) = draw_deck;
end

cluster_matrix_temp = cluster_matrix_L;

index_temp_1 = cluster_matrix_temp(:,1)./(cluster_matrix_temp(:,1) + cluster_matrix_temp(:,2));
index_temp_1(isnan(index_temp_1)) = 0.5;
index_temp_2 = cluster_matrix_temp(:,1) - cluster_matrix_temp(:,2);
index_temp = [index_temp_1, index_temp_2];

[~,index_order] = sortrows(index_temp,[1,2]);

cluster_lat = cluster_lat+0.5;
for draw_deck = 1:max(cluster_lat)-0.5
    cluster_lat(cluster_lat==index_order(draw_deck)+0.5) = draw_deck;
end

%%
cluster_geo_temp = cluster_geo;
cluster_dyn_temp = cluster_dyn;
scatter_size = 0.28;
mult_c = 2.5;

cluster_matrix_GD = zeros(max(cluster_geo_temp),max(cluster_dyn_temp),2);
for ss = 1:length(train_ROI)
    temp_type = (train_ROI(ss)>boundary)+1;
    cluster_matrix_GD(cluster_geo_temp(train_ROI(ss)), cluster_dyn_temp(train_ROI(ss)), temp_type) =  cluster_matrix_GD(cluster_geo_temp(train_ROI(ss)), cluster_dyn_temp(train_ROI(ss)), temp_type)+1;
end
temp = sum(cluster_matrix_GD,3);
num_temp = max(temp(:))-1;
figure
hold on
for gg = 1:max(cluster_geo_temp)
    for dd = 1:max(cluster_dyn_temp)
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
hold on
% for ii = 1:max(cluster_geo)
%     plot([ii,ii],[0, max(cluster_dyn)],'color',[0.5 0.5 0.5],'linewidth',0.25);
% end
% for ii = 1:max(cluster_dyn)
%     plot([0, max(cluster_geo)],[ii,ii],'color',[0.5 0.5 0.5],'linewidth',0.25);
% end
xlim([0 max(cluster_geo_temp)])
ylim([0 max(cluster_dyn_temp)])
axis xy image off

%%
cluster_geo_temp = cluster_geo;
cluster_dyn_temp = cluster_dyn;
scatter_size = 0.28;
mult_c = 2.5;
cluster_matrix_GD = zeros(max(cluster_geo_temp),max(cluster_dyn_temp),2);
for ss = 1:length(test_ROI)
    temp_type = (test_ROI(ss)>boundary)+1;
    cluster_matrix_GD(cluster_geo_temp(test_ROI(ss)), cluster_dyn_temp(test_ROI(ss)), temp_type) =  cluster_matrix_GD(cluster_geo_temp(test_ROI(ss)), cluster_dyn_temp(test_ROI(ss)), temp_type)+1;
end
temp = sum(cluster_matrix_GD,3);
num_temp = max(temp(:))-1;
figure
hold on
for gg = 1:max(cluster_geo_temp)
    for dd = 1:max(cluster_dyn_temp)
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
hold on
% for ii = 1:max(cluster_geo)
%     plot([ii,ii],[0, max(cluster_dyn)],'color',[0.5 0.5 0.5],'linewidth',0.25);
% end
% for ii = 1:max(cluster_dyn)
%     plot([0, max(cluster_geo)],[ii,ii],'color',[0.5 0.5 0.5],'linewidth',0.25);
% end
xlim([0 max(cluster_geo_temp)])
ylim([0 max(cluster_dyn_temp)])
axis xy image off
