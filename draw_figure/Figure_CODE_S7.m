clear;
close all;
load('Rev_fig_4_cluster_oldWT4AD3_180409.mat')

%%
close all
input_error = lat_error_ROI;
cluster_gain = cluster_gain_lat;


Z = linkage(input_error,'average');
norm_c = -floor(min(log10(Z(:,3))));
Z(:,3) = log10(Z(:,3))+norm_c;
cutoff_list = Z(:,3);
[~,max_index] = max(cluster_gain);
cutoff = (cutoff_list(max_index)+cutoff_list(max_index+1))/2;


figure
[~,~,outperm] = dendrogram(Z,0,'colorthreshold',cutoff);
hold on
plot([0 4149],[cutoff cutoff],'r')
ylim([ceil(min(cutoff_list)) ceil(max(cutoff_list))])
set(gca, 'xtick',[0 length(cluster_geo)], 'ytick', [ceil(min(Z(:,3))), ceil(max(Z(:,3)))])
set(gca, 'XTick', [])
set(gca, 'YTick', [])

tick_temp = [];
for ii = ceil(min(cutoff_list)):ceil(max(cutoff_list))
    temp = log10(linspace(10^(ii),10^(ii+1),11));
    temp(end) = [];
    tick_temp = [tick_temp,temp];
end

h=figure;
plot(cutoff_list,cluster_gain)
hold on
plot([cutoff,cutoff],[floor(min(cluster_gain)) ceil(max(cluster_gain))],'r')
xlim([ceil(min(cutoff_list)) ceil(max(cutoff_list))])
set(gca, 'XTick', tick_temp)
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])

temp_input_error = input_error;
for ii = 1:size(input_error,1)
    temp_input_error(ii,ii) = nan;
end
figure
imagesc(log2(temp_input_error(outperm,outperm)))
temp = colormap();
colormap(flipud(temp));
axis xy image off
cluster_input = cluster(Z,'cutoff',cutoff,'criterion','distance');
hold on
cluster_temp = cluster_input(outperm);
for cc = 1:max(cluster_input)
    temp = find(cluster_temp==cc);
    stat_point = min(temp);
    end_point = max(temp);
    plot([stat_point end_point end_point stat_point stat_point]-0.5,...
        [stat_point stat_point end_point end_point stat_point]-0.5,'r','linewidth',1); 
end
set(gca, 'XTick', [])
set(gca, 'YTick', [])

    
    
    
    
    