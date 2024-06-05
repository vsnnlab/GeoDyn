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