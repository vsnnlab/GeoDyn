clear;
close all;
%%
ACT1 = zeros(101,101,50);
ACT2 = zeros(101,101,50);
ACT3 = zeros(101,101,50);
ACT4 = zeros(101,101,50);
ACT5 = zeros(101,101,50);
noise = 0;
%% ACT 1
[xx,yy] = meshgrid(-50:50,-50:50);
rr = abs(xx+yy*1i);
for tt = 1:50
    Amp = 2*sin(pi/100*tt) + 4;
    sig = 0*sin(pi/50*tt) + 15;
    temp = max(-(Amp/sig)*rr + Amp,0);
    size_lim = 15+sin(pi/50*tt)*0;
    %     temp = temp.*(rr<size_lim);
    temp = min(temp,4);
    ACT1(:,:,tt) = temp + rand(101,101)*noise-noise/2;
end

%% ACT 2
[xx,yy] = meshgrid(-50:50,-50:50);
rr = abs(xx+yy*1i);
for tt = 1:50
    %%% Same geometric profile
    r2 = 15+tt*0.2;
    r1 = sqrt(r2.^2 - 225);
    r0 = (r1+r2)/2;
    sig = (r2-r1)/2;
    
    Amp = 0*sin(pi/100*tt) + 4;
    temp = max(-(Amp/sig)*abs(rr-r0) + Amp,0);
    temp = temp/sum(temp(:))*sum(sum(ACT1(:,:,tt)));
    ACT2(:,:,tt) = temp + rand(101,101)*noise-noise/2;
end

%% ACT 3
[xx,yy] = meshgrid(-50:50,-50:50);
rr = abs(xx+yy*1i);
for tt = 1:50
    r2 = 15+tt*0.2;
    r1 = sqrt(r2.^2 - 225);
    r0 = (r1+r2)/2;
    Amp = 2*sin(pi/100*tt) + 4;
    temp = (15 - sqrt(4*r0.*abs(rr-r0)))*Amp/15;
    temp = min(temp,4);
    ACT3(:,:,tt) = temp + rand(101,101)*noise-noise/2;
end

%% ACT 4
[xx,yy] = meshgrid(-50:50,-50:50);
rr = abs(xx+yy*1i);
for tt = 1:50
    r2 = 15+0.1*tt;
    r1 = sqrt(r2.^2 - 225);
    r0 = (r1+r2)/2;
    Amp = 2*sin(pi/100*tt) + 4;
    temp = (15 - sqrt(4*r0.*abs(rr-r0)))*Amp/15;
    temp = min(temp,4);
    ACT4(:,:,tt) = temp + rand(101,101)*noise-noise/2;
end
%% ACT 5
[xx,yy] = meshgrid(-50:50,-50:50);
rr = abs(xx+yy*1i);
for tt = 1:50
    r1 = 8*sin(pi/100*tt);
    r2 = sqrt(r1.^2 + 225);
    r0 = (r1+r2)/2;
    Amp = 2*sin(pi/100*tt) + 4;
    temp = (15 - sqrt(4*r0.*abs(rr-r0)))*Amp/15;
    temp = min(temp,4);
    ACT5(:,:,tt) = temp + rand(101,101)*noise-noise/2;
end

%%
saveVideo = 0;
if saveVideo
    writerObj = VideoWriter(['Sample_Video/ari_amp.avi']);
    open(writerObj);
end
FIG = figure;
set(gcf,'position',[300 420 900 360])
for tt = 1:50
    drawnow;
    subplot(1,2,1)
    imagesc(ACT1(:,:,tt))
    caxis([0 4])
    colorbar;
    colormap(jet)
    axis xy image
    subplot(1,2,2)
    imagesc(ACT5(:,:,tt))
    caxis([0 4])
    colorbar;
    colormap(jet)
    axis xy image
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
%%
figure
subplot(1,3,1)
imagesc(ACT1(:,:,5));
axis xy image off
caxis([0 4])
subplot(1,3,2)
imagesc(ACT1(:,:,25));
axis xy image off
caxis([0 4])
subplot(1,3,3)
imagesc(ACT1(:,:,45));
axis xy image off
colormap(jet)
caxis([0 4])
figure
subplot(1,3,1)
imagesc(ACT2(:,:,5));
axis xy image off
caxis([0 4])
subplot(1,3,2)
imagesc(ACT2(:,:,25));
axis xy image off
caxis([0 4])
subplot(1,3,3)
imagesc(ACT2(:,:,45));
axis xy image off
colormap(jet)
caxis([0 4])
figure
subplot(1,3,1)
imagesc(ACT3(:,:,5));
axis xy image off
caxis([0 4])
subplot(1,3,2)
imagesc(ACT3(:,:,25));
axis xy image off
caxis([0 4])
subplot(1,3,3)
imagesc(ACT3(:,:,45));
axis xy image off
colormap(jet)
caxis([0 4])
figure
subplot(1,3,1)
imagesc(ACT4(:,:,5));
axis xy image off
caxis([0 4])
subplot(1,3,2)
imagesc(ACT4(:,:,25));
axis xy image off
caxis([0 4])
subplot(1,3,3)
imagesc(ACT4(:,:,45));
axis xy image off
colormap(jet)
caxis([0 4])
%%
figure
subplot(1,3,1)
geo_temp1 = geo_profile(ACT1,ones(101),0);
imagesc(geo_temp1);
colormap(hot)
axis xy image
caxis([0 0.1])
subplot(1,3,2)
geo_temp2 = geo_profile(ACT4,ones(101),0);
imagesc(geo_temp2);
colormap(hot)
axis xy image
caxis([0 0.1])
hh= subplot(1,3,3);
imagesc(geo_temp2-geo_temp1)
colormap(hh,jet)
axis xy image
caxis([-0.02 0.02])
load('rbcmap.mat')

figure
subplot(1,3,1)
dyn_temp1 = dyn_profile(ACT1,ones(101));
imagesc(dyn_temp1);
colormap(hot)
axis xy image
caxis([-0.01 0.01])
subplot(1,3,2)
dyn_temp2 = dyn_profile(ACT4,ones(101));
imagesc(dyn_temp2);
colormap(rbcmap)
axis xy image
caxis([-max(dyn_temp2(:)) max(dyn_temp2(:))])
hh= subplot(1,3,3);
imagesc(dyn_temp2-dyn_temp1)
colormap(hh,jet)
axis xy image
caxis([-max(dyn_temp2(:)) max(dyn_temp2(:))])