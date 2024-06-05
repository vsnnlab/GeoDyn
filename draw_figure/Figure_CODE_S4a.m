close all;
clear;

%%
load('VSDI_prediction_total_WT6_AD4_011618.mat');
WTAD = prediction_rate_GD(:,end-2);
load('VSDI_prediction_AD1_AD4_011618.mat');
AD1AD4 = prediction_rate_GD(:,end-2);
load('VSDI_prediction_AD3_AD4_011618.mat');
AD3AD4 = prediction_rate_GD(:,end);
load('VSDI_prediction_WT1_WT5_011618.mat');
WT1WT5 = prediction_rate_GD(:,end);
load('VSDI_prediction_WT1_WT3_011618.mat');
WT1WT3 = prediction_rate_GD(:,end);

figure
hold on
bar([1 2 3 4 5],[mean(WTAD) mean(AD1AD4) mean(AD3AD4) mean(WT1WT5) mean(WT1WT3)],'facecolor',[0.75 0.75 0.75])
scatter(ones(100,1),WTAD);
scatter(ones(100,1)*2,AD1AD4);
scatter(ones(100,1)*3,AD3AD4);
scatter(ones(100,1)*4,WT1WT5);
scatter(ones(100,1)*5,WT1WT3);
xlim([0 6]);
ylim([0.5 0.9]);

figure
hold on
bar([1 2 3 4 5],[mean(WTAD) mean(AD1AD4) mean(AD3AD4) mean(WT1WT5) mean(WT1WT3)],'facecolor',[0.75 0.75 0.75])
errorbar([mean(WTAD) mean(AD1AD4) mean(AD3AD4) mean(WT1WT5) mean(WT1WT3)],[std(WTAD) std(AD1AD4) std(AD3AD4) std(WT1WT5) std(WT1WT3)]/10,'.')
xlim([0 6]);
ylim([0.5 0.8]);

figure
bar([mean(prediction_rate_GD(:,end)), mean(prediction_rate_ML(:,end)), mean(prediction_rate_M(:,end)), mean(prediction_rate_L(:,end))])
ylim([0.5 0.8])
hold on
errorbar([mean(prediction_rate_GD(:,end)), mean(prediction_rate_ML(:,end)), mean(prediction_rate_M(:,end)), mean(prediction_rate_L(:,end))]...
    ,[std(prediction_rate_GD(:,end)), std(prediction_rate_ML(:,end)), std(prediction_rate_M(:,end)), std(prediction_rate_L(:,end))]/10,'.')

