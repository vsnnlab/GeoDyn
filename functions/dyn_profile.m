function [P,P_norm] = dyn_profile(sample,mask);
% Smoothing parameters
alpha = 20;
beta = 0.01;

% define cordinates
[xx,yy] = meshgrid(1:size(sample,2),1:size(sample,1));
xx_temp = xx - (size(sample,2)+1)/2;
yy_temp = yy - (size(sample,1)+1)/2;
rr = abs(xx_temp + 1i*yy_temp);
rr_limit = min(size(sample,2),size(sample,1))/3;


%pre-define variable sizes
PVF = zeros(size(sample));
CM_xy = zeros(size(sample,1),2);

% Generate angular parameters
angle_dist = 7.5*pi/180;
angle_list = -pi:angle_dist/2:pi;
angle_vec = [cos(angle_list); sin(angle_list)];
angle_filt = zeros(size(sample,1),size(sample,2),length(angle_list));

P = zeros(length(angle_list),size(sample,3)-1);
P_norm = zeros(length(angle_list),size(sample,3)-1);

for tt = 1:size(sample,3)-1;
    frame1 = sample(:,:,tt).*mask;
    frame2 = sample(:,:,tt+1).*mask;
    
    % Calculate Center of Mass
    temp = frame1.*xx_temp;
    CM_xy(tt,1) = sum(temp(:))/sum(abs(frame1(:))) + (size(sample,2)+1)/2;
    temp = frame1.*yy_temp;
    CM_xy(tt,2) = sum(temp(:))/sum(abs(frame1(:))) + (size(sample,1)+1)/2;
    
    % Calculate Vector Field: Optic-Flow method
    [PVF_uu, PVF_vv] = PVF_search2(frame1, frame2, alpha,beta);
    PVF(:,:,tt) = PVF_uu + 1i*PVF_vv;
    
    %%% Absolute scale
    % Continueous gaussain angular mask
    
    for aa = 1:length(angle_list)
        theta = angle_norm(atan2(yy-CM_xy(tt,2),xx-CM_xy(tt,1))-angle_list(aa));
        angle_filt_temp = exp(-(theta.^2)/2/(angle_dist/2).^2).*(rr<rr_limit);
        angle_filt(:,:,aa) = angle_filt_temp;
    end
    
    weighted_PVF = repmat(PVF(:,:,tt),1,1,size(angle_filt,3)).*angle_filt;
    
    % normalize by average
    PVF_sum = squeeze(nansum(nansum(weighted_PVF,2),1));
    PVF_sum_vec = [real(PVF_sum)';imag(PVF_sum)'];
    projected_PVF_sum_vec = dot(PVF_sum_vec,angle_vec);
    projected_PVF_sum_vec(isnan(projected_PVF_sum_vec)) = 0;
    
    if tt == 1
        ref_angle = sum(angle_list.*abs(projected_PVF_sum_vec))/sum(abs(projected_PVF_sum_vec));
    end
    
    %%% normalized scale
    angle_list_norm = angle_norm(angle_list+ref_angle);
    angle_vec_norm = [cos(angle_list_norm); sin(angle_list_norm)];
    %     angle_list_norm
    % Continueous gaussain angular mask
    for aa = 1:length(angle_list_norm)
        theta = angle_norm(atan2(yy-CM_xy(tt,2),xx-CM_xy(tt,1))-angle_list_norm(aa));
        angle_filt_temp = exp(-(theta.^2)/2/(angle_dist/2).^2).*(rr<rr_limit);
        angle_filt(:,:,aa) = angle_filt_temp;
    end
    
    weighted_PVF = repmat(PVF(:,:,tt),1,1,size(angle_filt,3)).*angle_filt;
    
    % normalize by average
    PVF_sum = squeeze(nansum(nansum(weighted_PVF,2),1));
    PVF_sum_vec = [real(PVF_sum)';imag(PVF_sum)'];
    projected_PVF_sum_vec_norm = dot(PVF_sum_vec,angle_vec_norm);
    projected_PVF_sum_vec_norm(isnan(projected_PVF_sum_vec_norm)) = 0;
    
    P_norm(:,tt) = projected_PVF_sum_vec_norm;
    P(:,tt) = projected_PVF_sum_vec;
    
    ref_angle = sum(angle_list.*abs(projected_PVF_sum_vec))/sum(abs(projected_PVF_sum_vec));
    
end
end

%% Functions
function [uu, vv] = PVF_search2(im1, im2, alpha,beta, ite, uInitial, vInitial)
% Default parameters
if nargin<3
    alpha=1;
end
if nargin<4
    beta=1;
end
if nargin<5
    ite=100;
end
if nargin<6 || nargin<7
    uInitial = zeros(size(im1(:,:,1)));
    vInitial = zeros(size(im2(:,:,1)));
elseif size(uInitial,1) ==0 || size(vInitial,1)==0
    uInitial = zeros(size(im1(:,:,1)));
    vInitial = zeros(size(im2(:,:,1)));
end
% Convert images to grayscale
if size(size(im1),2)==3
    im1=rgb2gray(im1);
end
if size(size(im2),2)==3
    im2=rgb2gray(im2);
end
im1=double(im1);
im2=double(im2);

im1=smoothImg(im1,1);
im2=smoothImg(im2,1);
%
% Set initial value for the flow vectors
uu = uInitial;
vv = vInitial;

% Estimate spatiotemporal derivatives
[fx, fy, ft] = computeDerivatives(im1, im2);

% Kernels
kernel_1=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
jacob_ss_xx = [-1 1];
jacob_ss_yy = [-1; 1];

% Iterations
for i=1:ite
    % Compute local averages of the flow vectors
    uAvg=conv2(uu,kernel_1,'same');
    vAvg=conv2(vv,kernel_1,'same');
    
    %     error_p = ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft );
    error_s = convn(uu,jacob_ss_xx,'same').^2 + convn(uu,jacob_ss_yy,'same').^2 + convn(vv,jacob_ss_xx,'same').^2 + convn(vv,jacob_ss_yy,'same').^2;
    
    %     % Original: Compute flow vectors constrained by its local average and the optical flow constraints
    uu = uAvg - ( fx .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha *( beta^2 + error_s) );
    vv = vAvg - ( fy .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha *( beta^2 + error_s) );
    
    % Original: Compute flow vectors constrained by its local average and the optical flow constraints
    uu = uAvg - ( fx .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha^2 + fx.^2 + fy.^2);
    vv = vAvg - ( fy .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha^2 + fx.^2 + fy.^2);
end
uu(isnan(uu))=0;
vv(isnan(vv))=0;
end

function G=gaussFilter(segma,kSize)
if nargin<1
    segma=1;
end
if nargin<2
    kSize=2*(segma*3);
end

x=-(kSize/2):(1+1/kSize):(kSize/2);
G=(1/(sqrt(2*pi)*segma)) * exp (-(x.^2)/(2*segma^2));
end

function smoothedImg=smoothImg(img,segma)
if nargin<2
    segma=1;
end

G=gaussFilter(segma);
smoothedImg=conv2(img,G,'same');
smoothedImg=conv2(smoothedImg,G','same');
end


function [fx, fy, ft] = computeDerivatives(im1, im2)
if size(im2,1)==0
    im2=zeros(size(im1));
end
% Horn-Schunck original method
fx = conv2(im1,0.25* [-1 1; -1 1],'same') + conv2(im2, 0.25*[-1 1; -1 1],'same');
fy = conv2(im1, 0.25*[-1 -1; 1 1], 'same') + conv2(im2, 0.25*[-1 -1; 1 1], 'same');
ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
end



