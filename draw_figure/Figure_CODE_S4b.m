

cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
for ss = 1:length(ROI)
    if (ROI(ss)>421 && ROI(ss)<=1089)
        temp_type = 1;
        cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;

    elseif (ROI(ss)>1089 && ROI(ss)<=1548)
        temp_type = 2;
        cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
    end
end

%%%
scatter_size = 0.12;
mult_c = 5;
temp = sum(cluster_matrix_GD,3);
num_temp = max(temp(:))-1;
color_index = [ 0 255 0; 0 102 153; 255 180 0; 255 0 0]/255;
color_index = [255 180 0; 255 0 0]/255;
% color_index = flipud(color_index);
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
            for ii = 2:-1:1
                secdraw(gg,dd,0,360*cum_temp_ratio(ii),plot_size,color_index(ii,:));
            end
        end
    end
end

axis xy image
xlim([0 max(cluster_geo)+1])
ylim([0 max(cluster_dyn)+1])
%%
temp = squeeze(sum(cluster_matrix_GD,2));
figure
bar(temp)

%%
color_index = [ 0 255 0; 0 102 153; 255 180 0; 255 0 0]/255; 
%
cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),4);
for ss = 1:length(ROI)
    if (ROI(ss)>0 && ROI(ss)<=235)
        temp_type = 1;
    elseif (ROI(ss)>235 && ROI(ss)<=421)
        temp_type = 2;
    elseif (ROI(ss)>421 && ROI(ss)<=1089)
        temp_type = 3;
    elseif (ROI(ss)>1089 && ROI(ss)<=1548)
        temp_type = 4;
    end
    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
end

for individual = 1:4
    
    if individual >2
        indi_color = [linspace(1,0.9176,64)',linspace(1,0.3333,64)',linspace(1,0.0784,64)'];
    else
        indi_color = [linspace(1,0,64)',linspace(1,0.6,64)',linspace(1,0.2667,64)'];
    end
    
%     indi_color = [linspace(1,color_index(individual,1),64)',linspace(1,color_index(individual,2),64)',...
%         linspace(1,color_index(individual,3),64)'];
    %%%
    temp_CM_GD = cluster_matrix_GD(:,:,individual);
    
    P_indi = sum(temp_CM_GD(:))/sum(cluster_matrix_GD(:));
    P_GD = sum(cluster_matrix_GD,3)'/sum(cluster_matrix_GD(:));
    
    P_GD_indi = temp_CM_GD'/sum(temp_CM_GD(:));
    
    P_indi_GD = P_GD_indi.*P_indi./P_GD;
    P_indi_GD(isnan(P_indi_GD)) = 0;

    figure
    imagesc(P_indi_GD)
    colormap(indi_color)
    axis xy image off
%     caxis([0 max(P_GD_AD(:))])
%     colorbar;

    %%%
    scatter_size = 0.15;
    mult_c = 2;
    temp = sum(cluster_matrix_GD,3);
    num_temp = max(temp(:))-1;
    
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
%         plot([ii,ii],[0, max(cluster_dyn)],'color',[0.5 0.5 0.5],'linewidth',0.25);
%     end
%     for ii = 1:max(cluster_dyn)
%         plot([0, max(cluster_geo)],[ii,ii],'color',[0.5 0.5 0.5],'linewidth',0.25);
%     end
%     xlim([0 max(cluster_geo)])
%     ylim([0 max(cluster_dyn)])
%     axis xy image off
end
