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

noise_geo = 3;
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

noise_dyn =3;
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


%%
saveVideo = 1;

if saveVideo
    writerObj = VideoWriter(['Sample_Video/arti_dyn_prop.avi']);
    open(writerObj);
end

FIG = figure;
set(gcf,'position',[300 420 900 360])
for tt = 1:50
    drawnow;
    subplot(1,2,1)
    imagesc(ACT_dyn_con(:,:,tt))
    caxis([0 5])
    h = colorbar;
    h.Label.String = 'Amplitude (z)';
    h.Label.FontSize =  12;
    colormap(jet)
    axis xy image off
    title('Sample #4')
    
    subplot(1,2,2)
    imagesc(ACT_dyn_prop(:,:,tt))
    caxis([0 5])
    h = colorbar;
    h.Label.String = 'Amplitude (z)';
    h.Label.FontSize =  12;
    colormap(jet)
    axis xy image off
    title('Sample #6')
    if saveVideo
        pause(0.02);
        frame = getframe(FIG);
        writeVideo(writerObj,frame);
        writeVideo(writerObj,frame);
        writeVideo(writerObj,frame);
    end
end
if saveVideo
    close(writerObj);
    close(FIG);
end