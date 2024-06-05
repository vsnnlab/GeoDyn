clear
close all
%%
geo_con_mat = cell(50,1);
geo_same_mat = cell(50,1);
geo_amp_mat = cell(50,1);
geo_size_mat = cell(50,1);
dyn_con_mat = cell(50,1);
dyn_same_mat = cell(50,1);
dyn_speed_mat = cell(50,1);
dyn_prop_mat = cell(50,1);

geo_error_same = zeros(50,1);
geo_error_amp = zeros(50,1);
geo_error_size = zeros(50,1);
dyn_error_same = zeros(50,1);
dyn_error_speed = zeros(50,1);
dyn_error_prop = zeros(50,1);

for iter = 1:50
    
    ACT_geo_con = zeros(101,101,50);
    ACT_geo_same = zeros(101,101,50);
    ACT_geo_amp = zeros(101,101,50);
    ACT_geo_size = zeros(101,101,50);
    ACT_dyn_con = zeros(101,101,50);
    ACT_dyn_same = zeros(101,101,50);
    ACT_dyn_speed = zeros(101,101,50);
    ACT_dyn_prop = zeros(101,101,50);
    
    %
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    sig = 20;
    
    noise_geo = 0;
    for tt = 1:50
        x0 = 0;
        y0 = 0;
        rr = abs((xx-x0)+(yy-y0)*1i);
        
        Amp = 0*sin(pi/50*tt) + 3;
        temp = max(-(Amp/sig)*rr + Amp,0);
        temp(temp<0.5) = 0;
        ACT_geo_con(:,:,tt) = temp + rand(101,101)*noise_geo-noise_geo/2;
        ACT_geo_same(:,:,tt) = temp + rand(101,101)*noise_geo-noise_geo/2;
        
        Amp = 0.5*sin(pi/50*tt) + 3;
        temp = max(-(Amp/sig)*rr + Amp,0);
        temp(temp<0.5) = 0;
        ACT_geo_amp(:,:,tt) = temp + rand(101,101)*noise_geo-noise_geo/2;
        
        Amp = 0*sin(pi/50*tt) + 3;
        temp = max(-(Amp/sig)*rr + Amp,0);
        size_lim = 0.5-sin(pi/50*tt)*0.5;
        temp(temp<size_lim) = 0;
        ACT_geo_size(:,:,tt) = temp + rand(101,101)*noise_geo-noise_geo/2;
    end
    
    noise_dyn = 0;
    sig = 20;
    Amp = 3;
    for tt = 1:50
        x0 = -25 + tt;
        y0 = 0;
        [xx,yy] = meshgrid(-50:50,-50:50);
        rr = abs((xx-x0)+(yy-y0)*1i);
        temp = max(-(Amp/sig)*rr + Amp,0);
        temp(temp<0.5) = 0;
        ACT_dyn_con(:,:,tt) = temp + rand(101,101)*noise_dyn-noise_dyn/2;
        ACT_dyn_same(:,:,tt) = temp + rand(101,101)*noise_dyn-noise_dyn/2;
        
        x0 = -25 + tt;
        y0 = 0;
        y_lim = 50 - 0.1*tt;
        [xx,yy] = meshgrid(-50:50,linspace(-y_lim,y_lim,101));
        rr = abs((xx-x0)+(yy-y0)*1i);
        temp = max(-(Amp/sig)*rr + Amp,0);
        temp(temp<0.5) = 0;
        ACT_dyn_prop(:,:,tt) = temp + rand(101,101)*noise_dyn-noise_dyn/2;
    end
    x0 = -25;
    for tt = 1:50
        x0 = x0 + (1+(sin(pi/50*tt)-0.6237)*0.1);
        y0 = 0;
        [xx,yy] = meshgrid(-50:50,-50:50);
        rr = abs((xx-x0)+(yy-y0)*1i);
        temp = max(-(Amp/sig)*rr + Amp,0);
        temp(temp<0.5) = 0;
        ACT_dyn_speed(:,:,tt) = temp + rand(101,101)*noise_dyn-noise_dyn/2;
    end
    
    %%%
    sample = ACT_dyn_prop;
    figure
    subplot(1,3,1)
    imagesc(sample(:,:,1));
    axis xy image off
    caxis([0 3]);
    colormap(jet)
    subplot(1,3,2)
    imagesc(sample(:,:,25));
    axis xy image off
    caxis([0 3]);
    subplot(1,3,3)
    imagesc(sample(:,:,50));
    axis xy image off
    caxis([0 3]);
    
    %%%
    
    
    %%
    [xx,yy,zz] = meshgrid(-10:10,-10:10,-10:10);
    ST_filt = exp(-((xx).^2 + (yy).^2 + (zz).^2)/8);
    ST_filt = ST_filt/sum(ST_filt(:));
    
    ACT_geo_con = convn(ACT_geo_con,ST_filt,'same');
    ACT_geo_same = convn(ACT_geo_same,ST_filt,'same');
    ACT_geo_amp = convn(ACT_geo_amp,ST_filt,'same');
    ACT_geo_size = convn(ACT_geo_size,ST_filt,'same');
    ACT_dyn_con = convn(ACT_dyn_con,ST_filt,'same');
    ACT_dyn_same = convn(ACT_dyn_same,ST_filt,'same');
    ACT_dyn_speed = convn(ACT_dyn_speed,ST_filt,'same');
    ACT_dyn_prop = convn(ACT_dyn_prop,ST_filt,'same');
    
    geo_con = geo_profile(ACT_geo_con,ones(101,101),0);
    geo_con = geo_con(3:41,:);
    geo_same = geo_profile(ACT_geo_same,ones(101,101),0);
    geo_same = geo_same(3:41,:);
    geo_amp = geo_profile(ACT_geo_amp,ones(101,101),0);
    geo_amp = geo_amp(3:41,:);
    geo_size = geo_profile(ACT_geo_size,ones(101,101),0);
    geo_size = geo_size(3:41,:);
    
    dyn_con = dyn_profile(ACT_dyn_con,ones(101,101));
    dyn_con = dyn_con(:,5:end-4);
    dyn_same = dyn_profile(ACT_dyn_same,ones(101,101));
    dyn_same = dyn_same(:,5:end-4);
    dyn_speed = dyn_profile(ACT_dyn_speed,ones(101,101));
    dyn_speed = dyn_speed(:,5:end-4);
    dyn_prop = dyn_profile(ACT_dyn_prop,ones(101,101));
    dyn_prop = dyn_prop(:,5:end-4);
    
    %%
    geo_con_mat{iter} = geo_con;
    geo_same_mat{iter} = geo_same;
    geo_amp_mat{iter} = geo_amp;
    geo_size_mat{iter} = geo_size;
    dyn_con_mat{iter} = dyn_con;
    dyn_same_mat{iter} = dyn_same;
    dyn_speed_mat{iter} = dyn_speed;
    dyn_prop_mat{iter} = dyn_prop;
    
    temp = (geo_con-geo_same).^2;
    geo_error_same(iter) = sum(temp(:));
    temp = (geo_con-geo_amp).^2;
    geo_error_amp(iter) = sum(temp(:));
    temp = (geo_con-geo_size).^2;
    geo_error_size(iter) = sum(temp(:));
    
    temp = (dyn_con-dyn_same).^2;
    dyn_error_same(iter) = sum(temp(:));
    temp = (dyn_con-dyn_speed).^2;
    dyn_error_speed(iter) = sum(temp(:));
    temp = (dyn_con-dyn_prop).^2;
    dyn_error_prop(iter) = sum(temp(:));
end
figure
mean_error = [mean(geo_error_same), mean(geo_error_size), mean(geo_error_amp)];
std_error = [std(geo_error_same), std(geo_error_size) , std(geo_error_amp)];
errorbar(mean_error,std_error,'.');
xlim([0 4]);

figure
mean_error = [mean(dyn_error_same), mean(dyn_error_speed), mean(dyn_error_prop)];
std_error = [std(dyn_error_same), std(dyn_error_speed), std(dyn_error_prop)];
errorbar(mean_error,std_error,'.');
xlim([0 4]);

[p,h] = ttest(dyn_error_same, dyn_error_prop)


%%
total_geo = [geo_con_mat;geo_same_mat;geo_size_mat;geo_amp_mat];
total_dyn = [dyn_con_mat;dyn_same_mat;dyn_speed_mat;dyn_prop_mat;];


geo_error = zeros(200);
dyn_error = zeros(200);
for cc = 1:200
    for rr = cc:200
        geo_cc = total_geo{cc,1};
        geo_rr = total_geo{rr,1};
        
        temp = (geo_cc-geo_rr).^2;
        geo_error(cc,rr) = sum(temp(:));
        geo_error(rr,cc) =  geo_error(cc,rr);
        
        dyn_cc = total_dyn{cc,1};
        dyn_rr = total_dyn{rr,1};
        
        temp = (dyn_cc-dyn_rr).^2;
        dyn_error(cc,rr) = sum(temp(:));
        dyn_error(rr,cc) =  dyn_error(cc,rr);
    end
end

geo_error = geo_error/max(geo_error(:));
dyn_error = dyn_error/max(dyn_error(:));
figure
imagesc(1-geo_error)
imagesc(1-dyn_error)

save Rev_supp_fig_1_011118.mat

