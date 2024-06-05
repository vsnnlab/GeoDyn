function TIA_index = geo_profile(f_sample,mask, ifplot)
thre_intensity = 0:0.1:5;
TIA_index = zeros(length(thre_intensity),size(f_sample,3));
for tt = 1:size(f_sample,3);
    temp_frame = f_sample(:,:,tt).*mask;
    for ii = 1:length(thre_intensity);
        temp_thre_frame = temp_frame>thre_intensity(ii);
        TIA_index(ii,tt) = sum(temp_thre_frame(:))/sum(mask(:));
    end
end
% fft_TIA = abs(fftshift(fft2(TIA_index)));
if ifplot
    figure;
    imagesc(TIA_index);
    caxis([0 0.3])
    axis xy image;
    colormap(jet);
    colorbar;
end
