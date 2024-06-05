clear
close all

%%
iter_num = 100;
geo_mat = cell(9,iter_num);
dyn_mat = cell(9,iter_num);
for iter = 1:iter_num
    
    ACT_cell = cell(9,1);
    
    ACT1 = zeros(101,101,50);
    ACT2 = zeros(101,101,50);
    ACT3 = zeros(101,101,50);
    ACT4 = zeros(101,101,50);
    ACT5 = zeros(101,101,50);
    ACT6 = zeros(101,101,50);
    ACT7 = zeros(101,101,50);
    ACT8 = zeros(101,101,50);
    ACT9 = zeros(101,101,50);
    noise = 0;
    
    %% ACT 1
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        Amp = 0*sin(pi/50*tt) + 3;
        sig = 0*sin(pi/50*tt) + 15;
        temp = max(-(Amp/sig)*rr + Amp,0);
%         size_lim = 15+sin(pi/50*tt)*0;
%         temp = temp.*(rr<size_lim);
        ACT1(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    ACT_cell{1,1} = ACT1;
    %% ACT 2
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        Amp = 0*sin(pi/50*tt) + 3;
        sig = 0*sin(pi/50*tt) + 20;
        temp = max(-(Amp/sig)*rr + Amp,0);
        size_lim = 10+sin(pi/50*tt)*10;
        temp = temp.*(rr<size_lim);
        ACT2(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    ACT_cell{2,1} = ACT2;
    %% ACT 3
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        Amp = 2*sin(pi/50*tt) + 2;
        sig = 0*sin(pi/50*tt) + 20;
        temp = max(-(Amp/sig)*rr + Amp,0);
        size_lim = 15+sin(pi/50*tt)*0;
        temp = temp.*(rr<size_lim);
        ACT3(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    ACT_cell{3,1} = ACT3;
    %% ACT 4
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        Amp = 4*sin(pi/100*tt) + 2;
        sig = 0*sin(pi/50*tt) + 20;
        temp = max(-(Amp/sig)*rr + Amp,0);
        size_lim = 15+sin(pi/50*tt)*0;
        temp = temp.*(rr<size_lim);
        temp = min(temp,4);
        ACT4(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    ACT_cell{4,1} = ACT4;
    %% ACT 5
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        r0 = tt*0.6;
        Amp = 0*sin(pi/100*tt) + 3;
        sig = 0*sin(pi/50*tt) + 10;
        temp = max(-(Amp/sig)*abs(rr-r0) + Amp,0);
        ACT5(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    ACT_cell{5,1} = ACT5;
    %% ACT 6
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        x0 = -25 + tt;
        y0 = 0;
        Amp = 0*sin(pi/100*tt) + 3;
        sig = 0*sin(pi/50*tt) + 15;
        rr = abs((xx-x0)+(yy-y0)*1i);
        temp = max(-(Amp/sig)*rr + Amp,0);
%         size_lim = 15+sin(pi/50*tt)*0;
%         temp = temp.*(rr<size_lim);
        ACT6(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    ACT_cell{6,1} = ACT6;
    %% ACT 7
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        x0 = -25 + tt;
        y0 = sin(pi/12.5*tt)*7.5;
        Amp = 0*sin(pi/100*tt) + 3;
        sig = 0*sin(pi/50*tt) + 15;
        rr = abs((xx-x0)+(yy-y0)*1i);
        temp = max(-(Amp/sig)*rr + Amp,0);
%         size_lim = 15+sin(pi/50*tt)*0;
%         temp = temp.*(rr<size_lim);
        ACT7(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    ACT_cell{7,1} = ACT7;
    %% ACT 8
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        for aa = -30:60:30
            x0 = -25 + cos(pi/180*aa)*tt;
            y0 = sin(pi/180*aa)*tt;
            Amp = 0*sin(pi/100*tt) + 3;
            sig = 0*sin(pi/50*tt) + 20;
            rr = abs((xx-x0)+(yy-y0)*1i);
            temp = max(-(Amp/sig)*rr + Amp,0);
            size_lim = 15+sin(pi/50*tt)*0;
            temp = temp.*(rr<size_lim);
            ACT8(:,:,tt) = max(ACT8(:,:,tt),temp);
        end
        ACT8(:,:,tt) = ACT8(:,:,tt) + rand(101,101)*noise-noise/2;
    end
    ACT_cell{8,1} = ACT8;
    %% ACT 9
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    for tt = 1:50
        for aa = -30:5:30
            x0 = -25 + cos(pi/180*aa)*tt;
            y0 = sin(pi/180*aa)*tt;
            Amp = 0*sin(pi/100*tt) + 3;
            sig = 0*sin(pi/50*tt) + 20;
            rr = abs((xx-x0)+(yy-y0)*1i);
            temp = max(-(Amp/sig)*rr + Amp,0);
            size_lim = 15+sin(pi/50*tt)*0;
            temp = temp.*(rr<size_lim);
            ACT9(:,:,tt) = max(ACT9(:,:,tt),temp);
        end
        ACT9(:,:,tt) = ACT9(:,:,tt) + rand(101,101)*noise-noise/2;
    end
    ACT_cell{9,1} = ACT9;
    %% Filtering
    [xx,yy,zz] = meshgrid(-10:10,-10:10,-10:10);
    ST_filt = exp(-((xx).^2 + (yy).^2 + (zz).^2)/8);
    ST_filt = ST_filt/sum(ST_filt(:));
    ACT_filt_cell = cell(9,1);
    for ii = 1:9
        ACT_filt_cell{ii,1} = convn(ACT_cell{ii,1},ST_filt,'same');
    end
    
    for ii = 1:9
        geo_temp = geo_profile(ACT_filt_cell{ii,1},ones(101),0);
        geo_temp = geo_temp(3:41,:);
        geo_mat{ii,iter} = geo_temp;
    end
    for ii = 1:9
        dyn_temp = dyn_profile(ACT_filt_cell{ii,1},ones(101));
        dyn_temp = dyn_temp(:,5:end-4);
        dyn_mat{ii,iter} = dyn_temp;
    end
end
%%
geo_mat_serial = [geo_mat(1,:),geo_mat(2,:),geo_mat(3,:),geo_mat(4,:)...
    ,geo_mat(5,:),geo_mat(6,:),geo_mat(7,:),geo_mat(8,:),geo_mat(9,:)];
dyn_mat_serial = [dyn_mat(1,:),dyn_mat(2,:),dyn_mat(3,:),dyn_mat(4,:)...
    ,dyn_mat(5,:),dyn_mat(6,:),dyn_mat(7,:),dyn_mat(8,:),dyn_mat(9,:)];
geo_error = zeros(9*iter_num);
dyn_error = zeros(9*iter_num);
for cc = 1:9*iter_num
    for rr = cc:9*iter_num
        geo_cc = geo_mat_serial{1,cc};
        geo_rr = geo_mat_serial{1,rr};
        
        temp = (geo_cc-geo_rr).^2;
        geo_error(cc,rr) = sum(temp(:));
        geo_error(rr,cc) =  geo_error(cc,rr);
        
        dyn_cc = dyn_mat_serial{1,cc};
        dyn_rr = dyn_mat_serial{1,rr};
        
        temp = (dyn_cc-dyn_rr).^2;
        dyn_error(cc,rr) = sum(temp(:));
        dyn_error(rr,cc) =  dyn_error(cc,rr);
    end
end

%%
geo_profile_mat = zeros(39,50,9*iter_num);
dyn_profile_mat = zeros(37,41,9*iter_num);
for ii = 1:9*iter_num
    geo_profile_mat(:,:,ii) = geo_mat_serial{1,ii};
    dyn_profile_mat(:,:,ii) = dyn_mat_serial{1,ii};
end

%%
geo_error = geo_error/max(geo_error(:));
dyn_error = dyn_error/max(dyn_error(:));

geo_corr = 1 - geo_error;
dyn_corr = 1 - dyn_error;

Z_geo = linkage(geo_error,'average');
Z_dyn = linkage(dyn_error,'average');

figure
[~,~,outperm_geo] = dendrogram(Z_geo,0);
figure
[~,~,outperm_dyn] = dendrogram(Z_dyn,0);

[cluster_geo,cutoff_geo,cluster_gain_geo] = optimal_cluster3(geo_error,geo_profile_mat);
[cluster_dyn,cutoff_dyn,cluster_gain_dyn] = optimal_cluster3(dyn_error,dyn_profile_mat);

figure
imagesc((geo_corr(outperm_geo,outperm_geo)))
axis image off

figure
imagesc((dyn_corr(outperm_dyn,outperm_dyn)).^20)
axis image off

%%
close all
load('rbcmap.mat')
for iter = 1:9
    figure
%     h1 = subplot(1,2,1);
%     geo_avg = mean(geo_profile_mat(:,:,(iter-1)*100+1:(iter)*100),3);
%     imagesc(geo_avg);
%     colormap(h1,hot)
%     axis xy off;
%     caxis([0 0.15])
%     
%     h2 = subplot(1,2,2);
    dyn_avg = mean(dyn_profile_mat(:,:,(iter-1)*100+1:(iter)*100),3);
    imagesc(dyn_avg);
%     colormap(h2,rbcmap)
colormap(rbcmap)
    axis xy off;
    caxis([-0.3 0.3])
end



