%%%
cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),8);
for ss = 1:1000
    if (ROI(ss)>0 && ROI(ss)<=415)
        temp_type = 1;
    elseif (ROI(ss)>415 && ROI(ss)<=914)
        temp_type = 2;
    elseif (ROI(ss)>914 && ROI(ss)<=1567)
        temp_type = 3;
    elseif (ROI(ss)>1567 && ROI(ss)<=1959)
        temp_type = 4;
    elseif (ROI(ss)>1959 && ROI(ss)<=2244)
        temp_type = 5;
    elseif (ROI(ss)>2244 && ROI(ss)<=2687)
        temp_type = 6;
    elseif (ROI(ss)>2687 && ROI(ss)<=3327)
        temp_type = 7;
    else (ROI(ss)>3327 && ROI(ss)<=3897)
        temp_type = 8;
    end
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end

discrim_ratio = max(cluster_matrix_GD,[],3)./sum(cluster_matrix_GD,3);
discrim_ratio(isnan(discrim_ratio)) = 0;
figure
imagesc(discrim_ratio)
%%%
%%


scatter_size = 0.2;
mult_c = 2;
temp = sum(cluster_matrix_GD,3);
num_temp = max(temp(:))-1;
color_index = [ 0 255 0; 0 204 51; 0 153 102; 0 102 153; 0 51 204; 0 0 255; 255 145 0; 255 0 0]/255;
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
            for ii = 8:-1:1
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


%%
temp_index = find(ROI>3327 & ROI<=4000);
colorcode = color_index(8,:);
cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
for ss = temp_index
    temp_type = (ss>boundary)+1;
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end
%%%
P_AD = (length(temp_index)-boundary_temp)/length(temp_index);
P_WT = boundary_temp/length(temp_index);
P_GD = sum(cluster_matrix_GD,3)'/length(temp_index);
temp = cluster_matrix_GD(:,:,2)';
P_GD_AD = temp/sum(temp(:));
temp = cluster_matrix_GD(:,:,1)';
P_GD_WT = temp/sum(temp(:));
P_AD_GD = P_GD_AD.*P_AD./P_GD;
P_WT_GD = P_GD_WT.*P_WT./P_GD;

P_WTAD = (P_AD_GD-P_WT_GD)./(P_AD_GD+P_WT_GD);
P_WTAD(isnan(P_WTAD)) = 0;


%%%
scatter_size = 0.15;
mult_c = 2;
temp = sum(cluster_matrix_GD,3);
num_temp = max(temp(:))-1;
num_temp = 243;

figure
hold on
for gg = 1:max(cluster_geo)
    for dd = 1:max(cluster_dyn)
        sample_size = sum(cluster_matrix_GD(gg,dd,:));
        plot_size = ((mult_c-1)/num_temp*sample_size + (num_temp-mult_c+1)/num_temp)*scatter_size;
        if sample_size ~= 0
            if cluster_matrix_GD(gg,dd,1)>cluster_matrix_GD(gg,dd,1)
                secdraw(gg-0.5,dd-0.5,90,360,plot_size,colorcode)
                temp_ratio = cluster_matrix_GD(gg,dd,2)/sample_size;
                if temp_ratio ~= 0
                    secdraw(gg-0.5,dd-0.5,90,360*temp_ratio,plot_size,colorcode)
                end
            else
                secdraw(gg-0.5,dd-0.5,90,360,plot_size,colorcode)
                temp_ratio = cluster_matrix_GD(gg,dd,1)/sample_size;
                if temp_ratio ~= 0
                    secdraw(gg-0.5,dd-0.5,90,360*temp_ratio,plot_size,colorcode)
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



