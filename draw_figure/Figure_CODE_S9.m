clear
close all

geo_WT_mat = zeros(39,50,50);
geo_AD_mat = zeros(39,50,50);
geo_con1_mat = zeros(39,50,50);
geo_con2_mat = zeros(39,50,50);

dyn_WT_mat = zeros(37,41,50);
dyn_AD_mat = zeros(37,41,50);
dyn_con1_mat = zeros(37,41,50);
dyn_con2_mat = zeros(37,41,50);

max_WT_mat = zeros(101,101,50);
max_AD_mat = zeros(101,101,50);
max_con1_mat = zeros(101,101,50);
max_con2_mat = zeros(101,101,50);

lat_WT_mat = zeros(101,101,50);
lat_AD_mat = zeros(101,101,50);
lat_con1_mat = zeros(101,101,50);
lat_con2_mat = zeros(101,101,50);

for iter = 1:50
    ACT_WT = zeros(101,101,50);
    ACT_AD = zeros(101,101,50);
    ACT_con_1 = zeros(101,101,50);
    ACT_con_2 = zeros(101,101,50);
    
    [xx,yy] = meshgrid(-50:50,-50:50);
    rr = abs(xx+yy*1i);
    noise = 3;
    for tt = 1:50
        Amp = 2*sin(pi/50*tt) + 1;
        sig = 15*sin(pi/50*tt) + 15;
        temp = max(-(Amp/sig)*rr + Amp,0);
        size_lim = 0.5-sin(pi/50*tt)*0;
        temp(temp<size_lim) = 0;
        ACT_WT(:,:,tt) = temp + rand(101,101)*noise-noise/2;
        
        Amp = 0.5*sin(pi/50*tt) + 0.5;
        sig = 15*sin(pi/50*tt) + 15;
        temp = max(-(Amp/sig)*rr + Amp,0);
        size_lim = 0-sin(pi/50*tt)*0;
        temp(temp<size_lim) = 0;
        ACT_AD(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    
    for tt = 1:50
        x0 = -25 + tt;
        y0 = 0;
        rr = abs((xx-x0)+(yy-y0)*1i);
        Amp = 0*sin(pi/50*tt) + 3;
        sig = 0*sin(pi/50*tt) + 20;
        temp = max(-(Amp/sig)*rr + Amp,0);
        size_lim = 0.5-sin(pi/50*tt)*0;
        temp(temp<size_lim) = 0;
        ACT_con_1(:,:,tt) = temp + rand(101,101)*noise-noise/2;
        
        x0 = 0;
        y0 = 0;
        rr = abs((xx-x0)+(yy-y0)*1i);
        
        r0 = 0+0.5*tt;
        Amp = 0*sin(pi/50*tt) + 3;
        sig = 0*sin(pi/50*tt) + 20;
        temp = max(-(Amp/sig)*abs(rr-r0) + Amp,0);
        size_lim = 0-sin(pi/50*tt)*0;
        temp(temp<size_lim) = 0;
        ACT_con_2(:,:,tt) = temp + rand(101,101)*noise-noise/2;
    end
    
    %%
    % saveVideo = 0;
    %
    % sample_1 = ACT_WT;
    % sample_2 = ACT_AD;
    % if saveVideo
    %     writerObj = VideoWriter(['Sample_Video/PLM_compare_1.avi']);
    %     open(writerObj);
    % end
    %
    % FIG = figure;
    % % set(gcf,'position',[300 420 900 360])
    % for tt = 1:50
    %     drawnow;
    %     subplot(1,2,1)
    %     imagesc(sample_1(:,:,tt))
    %     caxis([0 4])
    %     colorbar;
    %     colormap(jet)
    %     axis xy image
    %     title('Sample 1');
    %     subplot(1,2,2)
    %     imagesc(sample_2(:,:,tt))
    %     caxis([0 4])
    %     colorbar;
    %     colormap(jet)
    %     axis xy image
    %     title('Sample 2');
    %     if saveVideo
    %         pause(0.02);
    %         frame = getframe(FIG);
    %         writeVideo(writerObj,frame);
    %         writeVideo(writerObj,frame);
    %         writeVideo(writerObj,frame);
    %     end
    % end
    % if saveVideo
    %     close(writerObj);
    %     close(FIG);
    % end
    %
    % figure
    % subplot(1,3,1);
    % imagesc(sample_1(:,:,1))
    % colormap(jet)
    % axis xy image
    % caxis([0 4])
    % subplot(1,3,2);
    % imagesc(sample_1(:,:,25))
    % colormap(jet)
    % axis xy image
    % caxis([0 4])
    % subplot(1,3,3);
    % imagesc(sample_1(:,:,50))
    % colormap(jet)
    % axis xy image
    % caxis([0 4])
    %
    % figure
    % subplot(1,3,1);
    % imagesc(sample_2(:,:,1))
    % colormap(jet)
    % axis xy image
    % caxis([0 4])
    % subplot(1,3,2);
    % imagesc(sample_2(:,:,25))
    % colormap(jet)
    % caxis([0 4])
    % axis xy image
    % subplot(1,3,3);
    % imagesc(sample_2(:,:,50))
    % colormap(jet)
    % axis xy image
    % caxis([0 4])
    %%
    load('rbcmap.mat')
    
    [xx,yy,zz] = meshgrid(-5:5,-5:5,-5:5);
    ST_filt = exp(-((xx).^2 + (yy).^2 + (zz).^2)/2);
    ST_filt = ST_filt/sum(ST_filt(:));
    
    
    [xx,yy,zz] = meshgrid(-10:10,-10:10,-10:10);
    ST_filt = exp(-((xx).^2 + (yy).^2 + (zz).^2)/8);
    ST_filt = ST_filt/sum(ST_filt(:));
    
    ACT_WT = convn(ACT_WT,ST_filt,'same');
    ACT_AD = convn(ACT_AD,ST_filt,'same');
    ACT_con_1 = convn(ACT_con_1,ST_filt,'same');
    ACT_con_2 = convn(ACT_con_2,ST_filt,'same');
    
    
    geo_WT = geo_profile(ACT_WT,ones(101),0);
    geo_WT = geo_WT(3:41,:);
    geo_AD = geo_profile(ACT_AD,ones(101),0);
    geo_AD = geo_AD(3:41,:);
    dyn_WT = dyn_profile(ACT_WT,ones(101));
    dyn_WT = dyn_WT(:,5:end-4);
    dyn_AD = dyn_profile(ACT_AD,ones(101));
    dyn_AD = dyn_AD(:,5:end-4);
    lat_WT = latencymap2(ACT_WT);
    lat_AD = latencymap2(ACT_AD);
    max_WT = max(ACT_WT,[],3);
    max_AD = max(ACT_AD,[],3);
    
    geo_con1 = geo_profile(ACT_con_1,ones(101),0);
    geo_con1 = geo_con1(3:41,:);
    geo_con2 = geo_profile(ACT_con_2,ones(101),0);
    geo_con2 = geo_con2(3:41,:);
    dyn_con1 = dyn_profile(ACT_con_1,ones(101));
    dyn_con1 = dyn_con1(:,5:end-4);
    dyn_con2 = dyn_profile(ACT_con_2,ones(101));
    dyn_con2 = dyn_con2(:,5:end-4);
    lat_con1 = latencymap(ACT_con_1);
    lat_con2 = latencymap(ACT_con_2);
    max_con1 = max(ACT_con_1,[],3);
    max_con2 = max(ACT_con_2,[],3);
    
    geo_WT_mat(:,:,iter) = geo_WT;
    geo_AD_mat(:,:,iter) = geo_AD;
    geo_con1_mat(:,:,iter) = geo_con1;
    geo_con2_mat(:,:,iter) = geo_con2;
    
    dyn_WT_mat(:,:,iter) = dyn_WT;
    dyn_AD_mat(:,:,iter) = dyn_AD;
    dyn_con1_mat(:,:,iter) = dyn_con1;
    dyn_con2_mat(:,:,iter) = dyn_con2;
    
    max_WT_mat(:,:,iter) = max_WT;
    max_AD_mat(:,:,iter) = max_AD;
    max_con1_mat(:,:,iter) = max_con1;
    max_con2_mat(:,:,iter) = max_con2;
    
    lat_WT_mat(:,:,iter) = lat_WT;
    lat_AD_mat(:,:,iter) = lat_AD;
    lat_con1_mat(:,:,iter) = lat_con1;
    lat_con2_mat(:,:,iter) = lat_con2;
end

%%
geo_total = cat(3, geo_WT_mat,geo_AD_mat,geo_con1_mat,geo_con2_mat);
dyn_total = cat(3, dyn_WT_mat,dyn_AD_mat,dyn_con1_mat,dyn_con2_mat);
max_total = cat(3, max_WT_mat,max_AD_mat,max_con1_mat,max_con2_mat);
lat_total = cat(3, lat_WT_mat,lat_AD_mat,lat_con1_mat,lat_con2_mat);

geo_error = zeros(200);
dyn_error = zeros(200);
max_error = zeros(200);
lat_error = zeros(200);
for cc = 1:200
    for rr = cc:200
       	geo_cc = geo_total(:,:,cc);
        geo_rr = geo_total(:,:,rr);
        
        temp = (geo_cc-geo_rr).^2;
        geo_error(cc,rr) = sum(temp(:));
        geo_error(rr,cc) =  geo_error(cc,rr);
        
        dyn_cc = dyn_total(:,:,cc);
        dyn_rr = dyn_total(:,:,rr);
        
        temp = (dyn_cc-dyn_rr).^2;
        dyn_error(cc,rr) = sum(temp(:));
        dyn_error(rr,cc) =  dyn_error(cc,rr);
        
        max_cc = max_total(:,:,cc);
        max_rr = max_total(:,:,rr);
        
        temp = (max_cc-max_rr).^2;
        max_error(cc,rr) = sum(temp(:));
        max_error(rr,cc) =  max_error(cc,rr);
        
        lat_cc = lat_total(:,:,cc);
        lat_rr = lat_total(:,:,rr);
        
        temp = (lat_cc-lat_rr).^2;
        lat_error(cc,rr) = sum(temp(:));
        lat_error(rr,cc) =  lat_error(cc,rr);
    end
end
geo_error = geo_error/max(geo_error(:));
dyn_error = dyn_error/max(dyn_error(:));
max_error = max_error/max(max_error(:));
lat_error = lat_error/max(lat_error(:));

geo_cluster = optimal_cluster3(geo_error,geo_total);
figure
imagesc(1-geo_error(ROI,ROI))
axis image off;

figure
imagesc(1-dyn_error(ROI,ROI))
axis image off;

figure
imagesc(1-max_error(ROI,ROI))
axis image off;

figure
imagesc(1-lat_error(ROI,ROI))
axis image off;

temp = 53;
ROI(ROI==temp) = [];
%%
figure
h1 = subplot(2,2,1);
imagesc(geo_WT);
colormap(h1,hot);
caxis([0 0.3])
axis xy image
h2 = subplot(2,2,2);
imagesc(geo_AD);
colormap(h2,hot);
caxis([0 0.3])
axis xy image
h3 = subplot(2,2,3);
imagesc(dyn_WT);
colormap(h3,rbcmap);
caxis([-0.2 0.2])
axis xy image
h4 = subplot(2,2,4);
imagesc(dyn_AD);
colormap(h4,rbcmap);
caxis([-0.2 0.2])
axis xy image

figure
h1 = subplot(2,2,1);
imagesc(max(sample_1,[],3));
colormap(h1,jet);
caxis([0 4])
axis xy image
h2 = subplot(2,2,2);
imagesc(max(sample_2,[],3));
colormap(h2,jet);
caxis([0 4])
axis xy image
h3 = subplot(2,2,3);
imagesc(lat_WT);
colormap(h3,cool);
caxis([0 50])
axis xy image
h4 = subplot(2,2,4);
imagesc(lat_AD);
colormap(h4,cool);
caxis([0 50])
axis xy image

%%
phase_lat = zeros(101,101);
sample = ACT_con_2;
for cc = 1:101
    for rr = 1:101
        temp_signal = squeeze(sample(cc,rr,:));
        complex_signal = hilbert(temp_signal);
        Amp = abs(complex_signal);
        phase = angle(complex_signal);
        lat = find(phase == 0);
        lat(lat==1) = [];
        lat(lat==50) = [];
        if ~isempty(lat)
            if length(lat) ==1
                phase_lat(cc,rr) = lat;
            else
                phase_lat(cc,rr) = 0;
            end
        else
            sign_phase = sign(phase);
            lat = find(abs(diff(sign_phase)) ==2)+1;
            trans_point = find(abs(diff(phase))>2.5)+1;
            for ii = 1:length(trans_point)
                lat(lat == trans_point(ii)) =[];
            end
            lat(lat==1) = [];
            lat(lat==50) = [];
            if length(lat) ==1
                phase_lat(cc,rr) = lat;
            else
                phase_lat(cc,rr) = 0;
            end
        end
    end
end

figure
imagesc(phase_lat)
axis xy image
colormap(cool)
caxis([0 50])