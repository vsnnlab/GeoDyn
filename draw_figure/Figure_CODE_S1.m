[cluster_geo,cutoff_geo,cluster_gain_geo,cutoff_list_geo] = optimal_cluster3(geo_error,geo_profile_mat);
[cluster_dyn,cutoff_dyn,cluster_gain_dyn,cutoff_list_dyn] = optimal_cluster3(dyn_error,dyn_profile_mat);

%%
Z = linkage(geo_error,'average');
Z(:,3) = log10(Z(:,3))+3;
cutoff_list_geo = Z(:,3);
 [~,max_index] = max(cluster_gain_geo);
cutoff_geo = (cutoff_list_geo(max_index)+cutoff_list_geo(max_index+1))/2;

figure
dendrogram(Z,0)
hold on
plot([0 900],[cutoff_geo cutoff_geo],'r')
ylim([0 5])
set(gca, 'xtick',[0 length(cluster_geo)], 'ytick', [ceil(min(Z(:,3))), ceil(max(Z(:,3)))])


tick_temp = [];
for ii = 0:5
    temp = log10(linspace(10^(ii),10^(ii+1),11));
    temp(end) = [];
    tick_temp = [tick_temp,temp];
end

h=figure;
plot(cutoff_list_geo,cluster_gain_geo)
hold on
plot([cutoff_geo,cutoff_geo],[0 1200],'r')
xlim([ceil(min(Z(:,3))), ceil(max(Z(:,3)))])
set(gca, 'XTick', tick_temp)

%%
Z = linkage(dyn_error,'average');
Z(:,3) = log10(Z(:,3))+3;
cutoff_list_dyn = Z(:,3);
 [~,max_index] = max(cluster_gain_dyn);
cutoff_dyn = (cutoff_list_dyn(max_index)+cutoff_list_dyn(max_index+1))/2;

figure
dendrogram(Z,0)
hold on
plot([0 900],[cutoff_dyn cutoff_dyn],'r')
ylim([0 5])
set(gca, 'xtick',[0 length(cluster_dyn)], 'ytick', [floor(min(cutoff_list_dyn)), ceil(max(cutoff_list_dyn))])


tick_temp = [];
for ii = 0:5
    temp = log10(linspace(10^(ii),10^(ii+1),11));
    temp(end) = [];
    tick_temp = [tick_temp,temp];
end

h=figure;
plot(cutoff_list_dyn,cluster_gain_dyn)
hold on
plot([cutoff_geo,cutoff_geo],[0 6000],'r')

%%
temp_index = [1:5,101:105,201:205,301:305,401:405,501:505,601:605,701:705,801:805];
dyn_error_Temp = dyn_error(temp_index,temp_index);

Z = linkage(dyn_error_Temp,'average');
Z(:,3) = log10(Z(:,3))+3;
figure
[~,~,outperm] = dendrogram(Z,0);
figure
imagesc((1-dyn_error_Temp(outperm,outperm)).^10)
axis image off
