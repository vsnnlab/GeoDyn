clear;
close all;

%%
load('filename.mat');

%%
f_bandpass = designfilt('bandpassiir','FilterOrder',8, ...
    'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',6, ...
    'SampleRate',150);

[xx,yy,zz] = meshgrid(-10:10,-10:10,-10:10);
ST_filt = exp(-((xx).^2 + (yy).^2 + (zz).^2)/8);
ST_filt = ST_filt/sum(ST_filt(:));

%%
for ff = 1:size(FOI,1)
    tic
    load(['data_oldWT4_180405/',FOI{ff}]);
   
    mean_c = repmat(mean(Trans_data_mat,3),1,1,3000);
    std_c = repmat(std(Trans_data_mat,[],3),1,1,3000);
    z_sample = (Trans_data_mat-mean_c)./std_c;
    
    perm_sample = permute(z_sample,[3,1,2]);
    
    f_perm_sample = filtfilt(f_bandpass, perm_sample);
    f_sample = permute(f_perm_sample,[2,3,1]);
    f_sample = gather(convn(gpuArray(f_sample), ST_filt,'same'));
    f_sample = f_sample(:,:,1:end);
    save(['filename',FOI{ff}],'f_sample');
    toc;
end