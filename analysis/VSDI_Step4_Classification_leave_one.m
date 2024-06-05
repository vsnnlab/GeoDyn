clear
close all
load('VSDI_profiles_oldWT4AD3_180408.mat');
load('VSDI_profile_error_oldWT4AD3_180408.mat');
%%
WT_ind = [1 235;236 421;422 544;545 622;];
AD_ind = [623 1081; 1082 1415; 1416 2083;];

%%
boundary = 622;
chance_level = 0.5;
train_num_set = linspace(30, 240, 8);
test_num = 60;
trial_num = 25;

%%
num_data = size(sample_profile,1);
sample_profile_mat_geo = zeros(51,30,num_data);
sample_profile_mat_dyn = zeros(97,29,num_data);
sample_profile_mat_max = zeros(126,97,num_data);
sample_profile_mat_lat = zeros(126,97,num_data);

for pp = 1:num_data
    sample_profile_mat_geo(:,:,pp) = imresize(sample_profile{pp,3},[51,30]);
    sample_profile_mat_dyn(:,:,pp) = imresize(sample_profile{pp,5},[97,29]);
    sample_profile_mat_max(:,:,pp) = sample_profile{pp,6};
    sample_profile_mat_lat(:,:,pp) = sample_profile{pp,7};
end

%%
result_mat = cell(4,3,6);
for wt = 1:4
    for ad = 1:3
        prediction_rate_GD = zeros(trial_num,length(train_num_set));
        prediction_rate_ML = zeros(trial_num,length(train_num_set));
        prediction_rate_G = zeros(trial_num,length(train_num_set));
        prediction_rate_D = zeros(trial_num,length(train_num_set));
        prediction_rate_M = zeros(trial_num,length(train_num_set));
        prediction_rate_L = zeros(trial_num,length(train_num_set));
        for train = 1:length(train_num_set)
            train_temp_num = train_num_set(train);
            for trial = 1:trial_num
                % Sample training & test set
                WT_length = WT_ind(wt,2) - WT_ind(wt,1) + 1; 
                temp_ind = randperm(WT_length);
                mother_WT = WT_ind(wt,1):WT_ind(wt,2);
                test_WT = mother_WT(temp_ind(1:test_num));
                
                AD_length = AD_ind(ad,2) - AD_ind(ad,1) + 1; 
                temp_ind = randperm(AD_length);
                mother_AD = AD_ind(wt,1):AD_ind(wt,2);
                test_AD = mother_AD(temp_ind(1:test_num));
                
                test_ROI = [test_WT,test_AD];
                
                mother_WT = 1:boundary;
                mother_WT(WT_ind(wt,1):WT_ind(wt,2)) = [];
                temp_ind = randperm(length(mother_WT));
                train_WT = mother_WT(temp_ind(1:train_temp_num));
                
                mother_AD = 1:size(sample_profile,1);
                mother_AD(AD_ind(ad,1):AD_ind(ad,2)) = [];
                mother_AD(1:boundary) = [];
                temp_ind = randperm(length(mother_AD));
                train_AD = mother_AD(temp_ind(1:train_temp_num));
                
                train_ROI = [train_WT, train_AD];
                
                % Clustering
                sim_set = [train_ROI, test_ROI];
                cluster_geo = optimal_cluster3(geo_error(sim_set,sim_set),sample_profile_mat_geo(:,:,sim_set));
                cluster_dyn = optimal_cluster3(dyn_error(sim_set,sim_set),sample_profile_mat_dyn(:,:,sim_set));
                cluster_max = optimal_cluster3(max_error(sim_set,sim_set),sample_profile_mat_max(:,:,sim_set));
                cluster_lat = optimal_cluster3(lat_error(sim_set,sim_set),sample_profile_mat_lat(:,:,sim_set));
                
                cluster_matrix_GD = zeros(max(cluster_geo),max(cluster_dyn),2);
                cluster_matrix_ML = zeros(max(cluster_max),max(cluster_lat),2);
                
                for ss = 1:length(train_ROI)
                    temp_type = (sim_set(ss)>boundary)+1;
                    cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type) =  cluster_matrix_GD(cluster_geo(ss), cluster_dyn(ss), temp_type)+1;
                    cluster_matrix_ML(cluster_max(ss), cluster_lat(ss), temp_type) =  cluster_matrix_ML(cluster_max(ss), cluster_lat(ss), temp_type)+1;
                end
                
                % sorting GD
                cluster_matrix_temp = sum(cluster_matrix_GD,2);
                index_temp_1 = cluster_matrix_temp(:,:,1)./(cluster_matrix_temp(:,:,1) + cluster_matrix_temp(:,:,2));
                index_temp_1(isnan(index_temp_1)) = 0.5;
                index_temp_2 = cluster_matrix_temp(:,:,1) - cluster_matrix_temp(:,:,2);
                index_temp = [index_temp_1, index_temp_2];
                [~,index_order] = sortrows(index_temp,[1,2]);
                cluster_geo = cluster_geo+0.5;
                for draw_deck = 1:max(cluster_geo)-0.5
                    cluster_geo(cluster_geo==index_order(draw_deck)+0.5) = draw_deck;
                end
                cluster_matrix_temp = sum(cluster_matrix_GD,1);
                index_temp_1 = cluster_matrix_temp(:,:,1)./(cluster_matrix_temp(:,:,1) + cluster_matrix_temp(:,:,2));
                index_temp_1(isnan(index_temp_1)) = 0.5;
                index_temp_2 = cluster_matrix_temp(:,:,1) - cluster_matrix_temp(:,:,2);
                index_temp = [index_temp_1', index_temp_2'];
                [~,index_order] = sortrows(index_temp,[1,2]);
                cluster_dyn = cluster_dyn+0.5;
                for draw_deck = 1:max(cluster_dyn)-0.5
                    cluster_dyn(cluster_dyn==index_order(draw_deck)+0.5) = draw_deck;
                end
                % sorting ML
                cluster_matrix_temp = sum(cluster_matrix_ML,2);
                index_temp_1 = cluster_matrix_temp(:,:,1)./(cluster_matrix_temp(:,:,1) + cluster_matrix_temp(:,:,2));
                index_temp_1(isnan(index_temp_1)) = 0.5;
                index_temp_2 = cluster_matrix_temp(:,:,1) - cluster_matrix_temp(:,:,2);
                index_temp = [index_temp_1, index_temp_2];
                [~,index_order] = sortrows(index_temp,[1,2]);
                cluster_max = cluster_max+0.5;
                for draw_deck = 1:max(cluster_max)-0.5
                    cluster_max(cluster_max==index_order(draw_deck)+0.5) = draw_deck;
                end
                cluster_matrix_temp = sum(cluster_matrix_ML,1);
                index_temp_1 = cluster_matrix_temp(:,:,1)./(cluster_matrix_temp(:,:,1) + cluster_matrix_temp(:,:,2));
                index_temp_1(isnan(index_temp_1)) = 0.5;
                index_temp_2 = cluster_matrix_temp(:,:,1) - cluster_matrix_temp(:,:,2);
                index_temp = [index_temp_1', index_temp_2'];
                [~,index_order] = sortrows(index_temp,[1,2]);
                cluster_lat = cluster_lat+0.5;
                for draw_deck = 1:max(cluster_lat)-0.5
                    cluster_lat(cluster_lat==index_order(draw_deck)+0.5) = draw_deck;
                end
                
                % classify GeoDyn
                train_data = [cluster_geo(1:length(train_ROI)),cluster_dyn(1:length(train_ROI))];
                class_data = ((sim_set(1:length(train_ROI))>boundary)+1)';
                svmStruct = fitcsvm(train_data,class_data);
                class_predict = predict(svmStruct,[cluster_geo(length(train_ROI)+1:end),cluster_dyn(length(train_ROI)+1:end)]);
                correct_bool_GD = class_predict== ((sim_set(length(train_ROI)+1:end)>boundary)+1)';
                
                % classify MAM + PLM
                train_data = [cluster_max(1:length(train_ROI)),cluster_lat(1:length(train_ROI))];
                class_data = ((sim_set(1:length(train_ROI))>boundary)+1)';
                svmStruct = fitcsvm(train_data,class_data);
                class_predict = predict(svmStruct,[cluster_max(length(train_ROI)+1:end),cluster_lat(length(train_ROI)+1:end)]);
                correct_bool_ML = class_predict== ((sim_set(length(train_ROI)+1:end)>boundary)+1)';
                
                % classify Geo
                train_data = cluster_geo(1:length(train_ROI));
                class_data = ((sim_set(1:length(train_ROI))>boundary)+1)';
                svmStruct = fitcsvm(train_data,class_data);
                class_predict = predict(svmStruct,cluster_geo(length(train_ROI)+1:end));
                correct_bool_G = class_predict==((sim_set(length(train_ROI)+1:end)>boundary)+1)';
                
                % classify Dyn
                train_data = cluster_dyn(1:length(train_ROI));
                class_data = ((sim_set(1:length(train_ROI))>boundary)+1)';
                svmStruct = fitcsvm(train_data,class_data);
                class_predict = predict(svmStruct,cluster_dyn(length(train_ROI)+1:end));
                correct_bool_D = class_predict==((sim_set(length(train_ROI)+1:end)>boundary)+1)';
                
                % classify MAM
                train_data = cluster_max(1:length(train_ROI));
                class_data = ((sim_set(1:length(train_ROI))>boundary)+1)';
                svmStruct = fitcsvm(train_data,class_data);
                class_predict = predict(svmStruct,cluster_max(length(train_ROI)+1:end));
                correct_bool_M = class_predict==((sim_set(length(train_ROI)+1:end)>boundary)+1)';
                
                % classify PLM
                train_data = cluster_lat(1:length(train_ROI));
                class_data = ((sim_set(1:length(train_ROI))>boundary)+1)';
                svmStruct = fitcsvm(train_data,class_data);
                class_predict = predict(svmStruct,cluster_lat(length(train_ROI)+1:end));
                correct_bool_L = class_predict==((sim_set(length(train_ROI)+1:end)>boundary)+1)';
                
                prediction_rate_GD(trial,train) = mean(correct_bool_GD);
                prediction_rate_ML(trial,train) = mean(correct_bool_ML);
                prediction_rate_G(trial,train) = mean(correct_bool_G);
                prediction_rate_D(trial,train) = mean(correct_bool_D);
                prediction_rate_M(trial,train) = mean(correct_bool_M);
                prediction_rate_L(trial,train) = mean(correct_bool_L);
            end
        end
        result_mat{wt,ad,1} = prediction_rate_GD;
        result_mat{wt,ad,2} = prediction_rate_ML;
        result_mat{wt,ad,3} = prediction_rate_G;
        result_mat{wt,ad,4} = prediction_rate_D;
        result_mat{wt,ad,5} = prediction_rate_M;
        result_mat{wt,ad,6} = prediction_rate_L;
    end
end

save VSDI_prediction_oldWT4AD3_040518.mat result_mat
    
%%
figure
count = 1;

prediction_rate_GD_mean = zeros(12,length(train_num_set));
prediction_rate_ML_mean = zeros(12,length(train_num_set));
prediction_rate_G_mean = zeros(12,length(train_num_set));
prediction_rate_D_mean = zeros(12,length(train_num_set));
prediction_rate_M_mean = zeros(12,length(train_num_set));
prediction_rate_L_mean = zeros(12,length(train_num_set));

prediction_rate_GD_std = zeros(12,length(train_num_set));
prediction_rate_ML_std = zeros(12,length(train_num_set));
prediction_rate_G_std = zeros(12,length(train_num_set));
prediction_rate_D_std = zeros(12,length(train_num_set));
prediction_rate_M_std = zeros(12,length(train_num_set));
prediction_rate_L_std = zeros(12,length(train_num_set));
        
for wt = 1:4
    for ad = 1:3
        subplot(3,4,count)
        prediction_rate_GD= result_mat{wt,ad,1};
        prediction_rate_ML = result_mat{wt,ad,2};
        prediction_rate_G = result_mat{wt,ad,3};
        prediction_rate_D= result_mat{wt,ad,4};
        prediction_rate_M = result_mat{wt,ad,5};
        prediction_rate_L = result_mat{wt,ad,6};
        
        for train = 1:length(train_num_set)
            temp = prediction_rate_GD(:,train,:);
            prediction_rate_GD_std(count,train) = std(temp(:))/sqrt(length(temp(:)));
            temp = prediction_rate_ML(:,train,:);
            prediction_rate_ML_std(count,train) = std(temp(:))/sqrt(length(temp(:)));
            temp = prediction_rate_G(:,train,:);
            prediction_rate_D_std(count,train) = std(temp(:))/sqrt(length(temp(:)));
            temp = prediction_rate_D(:,train,:);
            prediction_rate_D_std(count,train) = std(temp(:))/sqrt(length(temp(:)));
            temp = prediction_rate_M(:,train,:);
            prediction_rate_M_std(count,train) = std(temp(:))/sqrt(length(temp(:)));
            temp = prediction_rate_L(:,train,:);
            prediction_rate_L_std(count,train) = std(temp(:))/sqrt(length(temp(:)));
        end
        
        prediction_rate_GD_mean(count,:) = (mean(prediction_rate_GD,1));
        prediction_rate_ML_mean(count,:) = (mean(prediction_rate_ML,1));
        prediction_rate_G_mean(count,:) = (mean(prediction_rate_G,1));
        prediction_rate_D_mean(count,:) = (mean(prediction_rate_D,1));
        prediction_rate_M_mean(count,:) = (mean(prediction_rate_M,1));
        prediction_rate_L_mean(count,:) = (mean(prediction_rate_L,1));
        
        mean_pf = prediction_rate_GD_mean(count,:);
        std_pf = prediction_rate_GD_std(count,:);
        x = train_num_set;
        x = x';
        y = mean_pf';
        dy = std_pf';
        cc = [255 0 0]/255;
        
        figure; hold on;
        h=fill([x;flipud(x)],[y-dy;flipud(y+dy)],cc,'linestyle','none');
        set(h,'facealpha',0.2);
        plot(x,y,'.-','color',cc);
        
        mean_pf = prediction_rate_ML_mean(count,:);
        std_pf = prediction_rate_ML_std(count,:);
        y = mean_pf';
        dy = std_pf';
        cc = [220 220 0]/255;
        
        h=fill([x;flipud(x)],[y-dy;flipud(y+dy)],cc,'linestyle','none');
        set(h,'facealpha',0.2);
        plot(x,y,'.-','color',cc);
        
        count = count+1;
    end
end

