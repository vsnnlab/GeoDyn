function [min_error,ind,ourput_profile_1,output_profile_2] = profile_error(profile1,profile2)
%% Set buffer
buffer = 0;
%buffer = round(size(Prop_over,2)/4);

if (size(profile2,2)>=size(profile1,2))
    prop_ref = profile2;
    prop_over = profile1;
    order = 0;
else
    prop_ref = profile1;
    prop_over = profile2;
    order = 1;
end

%% Sliding on time axis
prop_ref_slide = zeros(size(prop_ref,1),size(prop_ref,2)+2*buffer);
prop_ref_slide(:,buffer+1:buffer+size(prop_ref,2)) = prop_ref;

mask_ref = zeros(size(prop_ref,1),size(prop_ref,2)+2*buffer);
mask_ref(:,buffer+1:buffer+size(prop_ref,2)) = 1;

error_array = zeros(1,size(prop_ref,2)-size(prop_over,2));
for slide_ii = 1:size(prop_ref,2)-size(prop_over,2)+1
    % padding profiles with zeros: matching size
    % give a shift with giving zeros
    prop_over_slide = zeros(size(prop_ref,1),size(prop_ref,2)+2*buffer);
    prop_over_slide(:,slide_ii:slide_ii+size(prop_over,2)-1) = prop_over;
    
    mask_over = zeros(size(prop_ref,1),size(prop_ref,2)+2*buffer);
    mask_over(:,slide_ii:slide_ii+size(prop_over,2)-1) = 1;
    % find the overlap region
    mask = (mask_ref.*mask_over)>0;
    % crop the non-overlap regions of the profiles
    prop_ref_crop = prop_ref_slide(mask);
    prop_ref_crop = reshape(prop_ref_crop,size(prop_ref,1),length(prop_ref_crop)/size(prop_ref,1));
    prop_over_crop = prop_over_slide(mask);
    prop_over_crop = reshape(prop_over_crop,size(prop_over,1),length(prop_over_crop)/size(prop_over,1));
    % measure the correlation
    temp = (prop_ref_crop-prop_over_crop).^2;
    error_array(slide_ii) = sqrt(sum(temp(:)));
end
[min_error,ind] = nanmin(error_array(:));

%% Generate Matched Profile

slide_ii = ind;
% padding profiles with zeros: matching size
% give a shift with giving zeros
prop_over_slide = zeros(size(prop_ref,1),size(prop_ref,2)+2*buffer);
prop_over_slide(:,slide_ii:slide_ii+size(prop_over,2)-1) = prop_over;
mask_over = zeros(size(prop_ref,1),size(prop_ref,2)+2*buffer);
mask_over(:,slide_ii:slide_ii+size(prop_over,2)-1) = 1;
% find the overlap region
mask = (mask_ref.*mask_over)>0;
% crop the non-overlap regions of the profiles
prop_ref_crop = prop_ref_slide(mask);
prop_ref_crop = reshape(prop_ref_crop,size(prop_ref,1),length(prop_ref_crop)/size(prop_ref,1));
prop_over_crop = prop_over_slide(mask);
prop_over_crop = reshape(prop_over_crop,size(prop_over,1),length(prop_over_crop)/size(prop_over,1));

if order
    ourput_profile_1 = prop_ref_crop;
    output_profile_2 = prop_over_crop;
else
    output_profile_2 = prop_ref_crop;
    ourput_profile_1 = prop_over_crop;
end
